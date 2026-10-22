## Run from anchored_v1: Rscript estimation/tests/check_firth_fixed_weights.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

num_jac <- function(f, z, h=1e-5) {
  sapply(seq_along(z), function(j) {
    e <- numeric(length(z)); e[j] <- h
    (f(z+e)-f(z-e))/(2*h)
  })
}
num_hessian_vector <- function(f, z, h=2e-4) {
  base <- f(z); k <- length(z)
  ans <- array(0, c(length(base), k, k))
  for (j in seq_len(k)) {
    ej <- numeric(k); ej[j] <- h
    ans[,j,j] <- (f(z+ej)-2*base+f(z-ej))/h^2
    if (j>1) for (l in seq_len(j-1)) {
      el <- numeric(k); el[l] <- h
      ans[,j,l] <- ans[,l,j] <-
        (f(z+ej+el)-f(z+ej-el)-f(z-ej+el)+f(z-ej-el))/(4*h^2)
    }
  }
  ans
}

## Independent stochastic-Taylor oracle using finite differences of the score.
## Expectations use the TRUE probabilities held fixed during differentiation.
## It does not use analytic H, C, or the simplified adjustment formula.
taylor_oracle <- function(theta, x, va, vb, w, deriv) {
  pa <- ncol(va); n <- length(w); k <- length(theta)
  aa <- seq_len(pa); bb <- pa + seq_len(ncol(vb))
  original <- deriv(theta[aa],theta[bb],x,va,vb)
  ptrue <- original$p
  contrib <- function(z, y) {
    t <- deriv(z[aa],z[bb],x,va,vb)
    D <- cbind(va*t$d.a,vb*t$d.b)
    D * (w*(y-t$p)/(t$p*(1-t$p)))
  }
  mean_score <- function(z) colSums(contrib(z,ptrue))
  J <- -num_jac(mean_score,theta)
  Ji <- solve(J)
  s1 <- contrib(theta,rep(1,n)); s0 <- contrib(theta,rep(0,n))
  contrast <- s1-s0
  K <- crossprod(contrast,contrast*(ptrue*(1-ptrue)))
  V <- Ji%*%K%*%t(Ji)
  covariance_term <- numeric(k)
  for (i in seq_len(n)) {
    ds <- num_jac(function(z) (contrib(z,rep(1,n))-contrib(z,rep(0,n)))[i,],theta)
    covariance_term <- covariance_term + ptrue[i]*(1-ptrue[i])*as.vector(ds%*%Ji%*%contrast[i,])
  }
  T <- num_hessian_vector(mean_score,theta)
  second_term <- .5*vapply(seq_len(k),function(j)sum(T[j,,]*V),numeric(1))
  list(adjustment=-covariance_term-second_term, J=J, V=V)
}


## Independent finite differences of the score, rather than the implemented
## probability-Hessian contractions. Also check the unit-weight Firth limit.
set.seed(20260920)
n <- 32
x <- rep(0:1, length.out = n)
z <- matrix(runif(n * 3, -.7, .7), n, 3)
vb <- cbind(1, z)
w <- exp(.8*z[,1]-.6*z[,2]+.3*z[,3]); w <- w/mean(w)
errors <- numeric()
for (model in c("RR", "RD")) for (pa in 1:3) {
  va <- cbind(1, z[, seq_len(pa-1), drop=FALSE])
  theta <- c(.2, rep(.12, pa-1), -1, .3, -.2, .1)
  deriv <- if (model == "RR") firth_rr_derivatives else firth_rd_derivatives
  terms <- deriv(theta[seq_len(pa)], theta[-seq_len(pa)], x, va, vb)
  y <- as.integer(seq_len(n) %% 3 == 0)
  got <- firth_score_adjustment(y, va, vb, w, terms)
  oracle <- taylor_oracle(theta, x, va, vb, w, deriv)
  errors <- c(errors, max(abs(got$adjustment-oracle$adjustment)))
  stopifnot(tail(errors, 1) < 5e-5,
            max(abs(got$covariance-oracle$V)) < 1e-6)
  unit <- firth_score_adjustment(y, va, vb, rep(1,n), terms)
  Iinv <- solve(unit$information)
  old <- numeric(length(theta))
  for (i in seq_len(n)) {
    D <- c(va[i,]*terms$d.a[i],vb[i,]*terms$d.b[i])
    H <- rbind(cbind(tcrossprod(va[i,])*terms$d.aa[i],
                      outer(va[i,],vb[i,])*terms$d.ab[i]),
               cbind(outer(vb[i,],va[i,])*terms$d.ab[i],
                      tcrossprod(vb[i,])*terms$d.bb[i]))
    old <- old + .5*D/(terms$p[i]*(1-terms$p[i]))*sum(diag(Iinv%*%H))
  }
  stopifnot(max(abs(unit$adjustment-old)) < 1e-9,
            max(abs(unit$weight.adjustment)) < 1e-9)
  for (mult in c(.01, 10, 1e4)) {
    scaled <- firth_score_adjustment(y, va, vb, mult*w, terms)
    stopifnot(max(abs(scaled$adjusted.score/mult-got$adjusted.score)) < 1e-9,
              max(abs(scaled$covariance-got$covariance)) < 1e-9)
  }
}
cat(sprintf("PASS: RR/RD stochastic-Taylor oracle (pa=1/2/3), maximum error %.3g\n", max(errors)))
cat("PASS: unit-weight Firth limit; score/covariance equivariance under weight rescaling\n")

## Independent two-arm delta-method oracle in the actual log RR / atanh RD
## and log OP parameterization. Each arm's weighted mean has known variance.
p0 <- .2; p1 <- .4
x <- rep(0:1, each=100)
w <- c(rep(c(.2,1.8),each=50), rep(c(.5,1.5),each=50))
va <- vb <- matrix(1,length(x),1)
qvar <- c(p0*(1-p0)*sum(w[x==0]^2)/sum(w[x==0])^2,
          p1*(1-p1)*sum(w[x==1]^2)/sum(w[x==1])^2)
for (model in c("RR", "RD")) {
  deriv <- if (model == "RR") firth_rr_derivatives else firth_rd_derivatives
  alpha <- if(model == "RR") log(p1/p0) else atanh(p1-p0)
  beta <- qlogis(p0)+qlogis(p1)
  terms <- deriv(alpha,beta,x,va,vb)
  got <- firth_score_adjustment(as.integer(seq_along(x)%%3==0),va,vb,w,terms)
  rho <- p1-p0
  bias.alpha <- if (model=="RR") .5*(qvar[1]/p0^2-qvar[2]/p1^2) else
    rho/(1-rho^2)^2*sum(qvar)
  bias.beta <- .5*sum(c((2*p0-1)/(p0*(1-p0))^2,
                        (2*p1-1)/(p1*(1-p1))^2)*qvar)
  stopifnot(max(abs(c(bias.alpha,bias.beta)+
                     solve(got$information,got$adjustment))) < 1e-10)
}
cat("PASS: two-arm RR/RD analytic leading-bias cancellation\n")

## Actual fits: nonuniform weights, RR/RD, with and without effect modifiers;
## verify final score independently of the solver's stored diagnostics.
dat <- expand.grid(id=1:100,x=0:1,z1=0:1)
rates <- c(.2,.3,.3,.6)
dat$y <- as.integer(dat$id <= 100*rates[1+dat$x+2*dat$z1])
ipd <- list(data=dat)
wf <- get.weight(covariate_matrix(dat),c(z1=.4))
for (model in c("RR","RD")) for (modifiers in list(NULL,"z1")) {
  args <- list(group="AD",data=ipd,weight=wf$weight,method="firth",
               effect_modifiers=modifiers,param=model)
  fit <- do.call(get.estimate,args)
  a <- seq_len(ncol(fit$va)); b <- length(a)+seq_len(ncol(fit$vb))
  deriv <- if(model=="RR") firth_rr_derivatives else firth_rd_derivatives
  terms <- deriv(fit$point.est[a],fit$point.est[b],dat$x,fit$va,fit$vb)
  check <- firth_score_adjustment(dat$y,fit$va,fit$vb,wf$weight,terms)
  stopifnot(isTRUE(fit$converged),
            max(abs(check$adjusted.score/fit$fc.diagnostics$score_scale)) <= 1e-6,
            max(abs(fit$cov-check$covariance)) < 1e-10,
            identical(fit$variance.method,"fixed_weight_model_sandwich_first_order"))
  args$weight <- 17*wf$weight
  scaled <- do.call(get.estimate,args)
  stopifnot(max(abs(fit$point.est-scaled$point.est)) < 1e-9,
            max(abs(fit$cov-scaled$cov)) < 1e-9)
  ## Covariance independently transformed from the arm-risk parameters for
  ## saturated modifier models (two distinct strata).
  if (length(modifiers)) {
    cell <- 1 + dat$x + 2*dat$z1
    risk <- vapply(1:4,function(k)terms$p[which(cell==k)[1]],numeric(1))
    risk.var <- vapply(1:4,function(k) {
      wk <- wf$weight[cell==k]
      risk[k]*(1-risk[k])*sum(wk^2)/sum(wk)^2
    },numeric(1))
    G <- matrix(0,4,4)
    d0 <- if(model=="RR") c(-1/risk[1],1/risk[2]) else
      c(-1,1)/(1-(risk[2]-risk[1])^2)
    d1 <- if(model=="RR") c(-1/risk[3],1/risk[4]) else
      c(-1,1)/(1-(risk[4]-risk[3])^2)
    G[1,1:2] <- d0; G[2,] <- c(-d0,d1)
    G[3,1:2] <- 1/(risk[1:2]*(1-risk[1:2]))
    G[4,] <- c(-G[3,1:2],1/(risk[3:4]*(1-risk[3:4])))
    oracle.cov <- G%*%diag(risk.var)%*%t(G)
    stopifnot(max(abs(fit$cov-oracle.cov)) < 1e-10)
  }
}
cat("PASS: RR/RD full fits, root residuals, fixed-weight covariance and scale invariance\n")

## Joint variance uses gamma matching the point-estimation weights.
for (model in c("RR","RD")) {
  args <- list(group="IPD",data=ipd,weight=wf$weight,gamma.ipd=wf$ga,
               target_mean=.4,method="firth",param=model)
  fit <- do.call(get.estimate,args)
  args$weight <- 13*wf$weight
  scaled <- do.call(get.estimate,args)
  stopifnot(max(abs(fit$point.est-scaled$point.est)) < 1e-9,
            max(abs(fit$cov-scaled$cov)) < 1e-9,
            identical(fit$variance.method,"joint_sandwich_first_order"))
  args$gamma.ipd <- 0
  bad <- tryCatch(do.call(get.estimate,args),error=identity)
  stopifnot(inherits(bad,"error"),grepl("matching the supplied weights",conditionMessage(bad)))
}
cat("PASS: joint sandwich scale invariance and mismatched-gamma rejection\n")

## Zero weights drop observations; a failed root must not become an available fit.
v <- cbind(1,dat$z1)
w <- rep(1,nrow(dat)); w[c(1,52,206,357)] <- 0
args <- list(param="RR",y=dat$y,x=dat$x,va=v,vb=v,weight=w,
             alpha.start=c(0,0),beta.start=c(0,0),pa=2,pb=2,max.step=100,thres=1e-6)
fit <- do.call(max.likelihood.firth.rr,args)
keep <- w > 0
small <- args
small$y <- dat$y[keep]; small$x <- dat$x[keep]
small$va <- small$vb <- v[keep,,drop=FALSE]; small$weight <- w[keep]
dropped <- do.call(max.likelihood.firth.rr,small)
stopifnot(isTRUE(fit$convergence),isTRUE(dropped$convergence),
          max(abs(fit$par-dropped$par)) < 1e-6)
args$max.step <- 1
bad <- do.call(max.likelihood.firth.rr,args)
stopifnot(!bad$convergence,bad$diagnostics$max_scaled_score > 1e-6)
bad <- tryCatch(MLEst(group="AD",param="RR",y=dat$y,x=dat$x,va=v,vb=v,
       weight=w,gamma.weight=0,max.step=1,thres=1e-6,alpha.start=c(0,0),
       beta.start=c(0,0),pa=2,pb=2,target_mean=0,method="firth"),error=identity)
stopifnot(isFALSE(bad$converged),all(is.finite(bad$point.est)),all(is.na(bad$se.est)),
          isFALSE(bad$last_optimizer$convergence),!is.null(bad$optimization_input))
args$max.step <- 100
args$vb <- cbind(v,v[,2]); args$pb <- 3; args$beta.start <- rep(0,3)
singular <- tryCatch(do.call(max.likelihood.firth.rr,args),error=identity)
stopifnot(inherits(singular,"error"))
cat("PASS: zero weights, failed-root rejection and singular-design rejection\n")
