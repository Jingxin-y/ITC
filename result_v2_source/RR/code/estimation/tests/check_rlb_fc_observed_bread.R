## Run from the project root. The oracle differentiates the ORIGINAL score,
## independently of glm_fc_score() and the analytic observed-bread expression.
source("load_anchored.R", encoding="UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R", encoding="UTF-8")
close_to <- function(a, b, tol=1e-7) stopifnot(
  identical(dim(a), dim(b)), max(abs(a-b)) < tol)
jacobian <- function(fun, a, step=1e-5) sapply(seq_along(a), function(j) {
  delta <- rep(0, length(a)); delta[j] <- step
  (fun(a+delta)-fun(a-delta))/(2*step)
})
dat <- expand.grid(id=1:50, x=0:1, z1=0:1, z2=0:1)
dat$y <- as.integer(dat$id <= 7+3*dat$x+4*dat$z1+2*dat$z2+2*dat$x*dat$z1)
cols <- c("z1", "z2"); Z <- as.matrix(dat[cols])
target <- c(z1=.4, z2=.6)
weights <- get.weight(Z, target)
v <- weights$weight/mean(weights$weight)
base <- list(data=dat, model_covariates=cols, effect_modifiers="z1", tolerance=1e-9)
for (w in list(rep(1,nrow(dat)), v)) {
  f <- suppressWarnings(do.call(fit_robust_lb_fc, c(base, list(w=w))))
  model <- suppressWarnings(do.call(fit_logbinomial_fc, c(base, list(w=w))))
  stopifnot(f$converged, model$converged, f$variance_status=="available",
    f$fit$variance_bread_type=="observed", model$fit$variance_bread_type=="expected")
  close_to(f$point.est, model$point.est, 1e-12)
  close_to(f$fc.diagnostics$adjusted_score, model$fc.diagnostics$adjusted_score, 1e-12)
  close_to(f$fit$information, model$fit$information, 1e-12)
  X <- model.matrix(f$fit); n <- nrow(X); theta <- coef(f$fit)
  score <- function(b) {
    p <- exp(drop(X%*%b))
    X*((dat$y-p)/(1-p))
  }
  H <- -jacobian(function(b) colMeans(w*score(b)), theta)
  Hi <- solve(H)
  oracle <- Hi%*%crossprod(w*score(theta))%*%t(Hi)/n^2
  treatment <- match("x", colnames(X))
  close_to(f$fit$variance_bread, H)
  close_to(f$vcov, oracle)
  close_to(f$se^2, oracle[treatment,treatment])
  close_to(c(f$low,f$up), f$est+c(-1,1)*qnorm(.975)*sqrt(oracle[treatment,treatment]))
  scaled <- suppressWarnings(do.call(fit_robust_lb_fc, c(base,list(w=17*w))))
  close_to(scaled$point.est,theta,1e-10); close_to(scaled$vcov,oracle)
  ## MAIC risk standardization with independently fixed model weights.
  d0 <- d1 <- dat; d0$x <- 0; d1$x <- 1
  X0 <- model.matrix(delete.response(terms(f$fit)),d0)
  X1 <- model.matrix(delete.response(terms(f$fit)),d1)
  target_fun <- function(b) log(weighted.mean(exp(drop(X1%*%b)),v))-
    log(weighted.mean(exp(drop(X0%*%b)),v))
  gradient <- drop(jacobian(target_fun, theta))
  terms <- maic_glm_terms(f$fit,dat,w)
  marginal <- maic_delta(terms,v,NULL,NULL,NULL,FALSE,w)
  close_to(marginal$se^2,drop(t(gradient)%*%oracle%*%gradient))
}
cat("PASS: unit/nonuniform weights, numerical observed-score derivative, coefficient and MAIC marginal Wald, FC point unchanged, weight scaling\n")

## Full stacked-score differentiation: calibration, regression, marginal target.
f <- suppressWarnings(do.call(fit_robust_lb_fc,c(base,list(w=v,
  gamma=weights$ga,target_mean=target,weight_covariates=cols,
  account_for_gamma=TRUE,account_for_target_mean=TRUE,
  target_mean_vcov=diag(c(.24,.24))/300))))
stopifnot(f$converged, f$variance_status=="available")
X <- model.matrix(f$fit); n <- nrow(X); k <- ncol(X); q <- ncol(Z)
theta <- coef(f$fit); D <- sweep(Z,2,target,"-")
d0 <- d1 <- dat; d0$x <- 0; d1$x <- 1
X0 <- model.matrix(delete.response(terms(f$fit)),d0)
X1 <- model.matrix(delete.response(terms(f$fit)),d1)
m0 <- weighted.mean(exp(drop(X0%*%theta)),v)
m1 <- weighted.mean(exp(drop(X1%*%theta)),v)
eta <- log(m1/m0)
stacked <- function(a) {
  w <- exp(drop(D%*%a[seq_len(q)])); w <- w/mean(w)
  b <- a[q+seq_len(k)]; p <- exp(drop(X%*%b))
  S <- X*((dat$y-p)/(1-p))
  final <- w*(exp(drop(X1%*%b))/m1-exp(drop(X0%*%b))/m0)-(a[q+k+1]-eta)
  cbind(w*D,w*S,final)
}
a <- c(weights$ga,theta,eta)
A <- jacobian(function(a) colMeans(stacked(a)),a); Ai <- solve(A)
influence <- -stacked(a)%*%t(Ai)
Vm <- diag(c(.24,.24))/300
ipd_cov <- crossprod(influence)/n^2
target_cov <- Ai[,seq_len(q),drop=FALSE]%*%Vm%*%t(Ai[,seq_len(q),drop=FALSE])
ix <- q+seq_len(k)
close_to(f$vcov.ipd,ipd_cov[ix,ix])
close_to(f$vcov.target,target_cov[ix,ix])
close_to(f$vcov,(ipd_cov+target_cov)[ix,ix])
terms <- maic_glm_terms(f$fit,dat,v)
marginal <- maic_delta(terms,v,Z,target,Vm,TRUE,v)
close_to(marginal$variance_ipd,ipd_cov[q+k+1,q+k+1])
close_to(marginal$variance_target,target_cov[q+k+1,q+k+1])
close_to(marginal$se^2,(ipd_cov+target_cov)[q+k+1,q+k+1])
old_fit <- f$fit; old_fit$variance_bread <- NULL
close_to(maic_glm_terms(old_fit,dat,v)$H,old_fit$information/n,1e-12)
cat("PASS: joint coefficient and marginal covariance match independently differentiated stacked equations, including AD means; legacy fit fallback\n")
