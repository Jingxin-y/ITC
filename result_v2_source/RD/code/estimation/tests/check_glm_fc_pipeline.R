source("load_anchored.R",encoding="UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R",encoding="UTF-8")
close_to <- function(a,b,tol=1e-6) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))
design <- webster_rr_design(event="rare",ess_ratio=.5)
spec <- make_brm_score_spec(); methods <- anchored_glm_fc_method_names()
for (mode in c("unit","maic")) {
  z <- suppressWarnings(run.anchored(design,500,seed=811,score_spec=spec,
    weight_covariates=paste0("z",1:6),model_weighting=mode,return_data=TRUE))
  stopifnot(setequal(names(z$fits),anchored_method_names()),length(z$fits)==18L,
    all(subset(z$estimates,method %in% methods)$status=="available"),
    all(subset(z$estimates,method %in% methods)$interval=="Wald"))
  for (j in seq(1,6,2)) close_to(z$fits[[methods[j]]]$point.est,z$fits[[methods[j+1]]]$point.est)
  dat <- apply_brm_scores(z$data$data.IPD,z$variables$score_spec)$data
  for (method in methods) {
    f <- z$fits[[method]]
    close_to(f$fit$prior.weights,z$model_weights$weight)
    if (mode=="unit") stopifnot(all(f$vcov.target==0))
  }
  ## Joint coefficient sandwich from a separately assembled stacked bread/meat.
  if (mode=="maic") {
    f <- z$fits[["RLP-BFC"]]; w <- z$model_weights$weight; n <- nrow(dat)
    X <- model.matrix(f$fit); p <- fitted(f$fit); S <- X*(dat$y-p)
    Z <- as.matrix(dat[paste0("z",1:6)]); mu <- z$data$data.AD$anchored.mean.x
    D <- sweep(Z,2,mu,"-"); Hg <- sweep(D,2,colSums(w*D)/sum(w),"-")
    Ag <- crossprod(D,w*Hg)/n; H <- crossprod(X,X*(w*p))/n; G <- crossprod(S,w*Hg)/n
    q <- ncol(Z); k <- ncol(X)
    A <- rbind(cbind(Ag,matrix(0,q,k)),cbind(G,-H))
    influence <- -(cbind(w*D,w*S))%*%t(MASS::ginv(A))
    sensitivity <- solve(H)%*%G%*%MASS::ginv(Ag)
    Vm <- z$data$data.AD$anchored.cov.x/z$data$data.AD$anchored.n
    expected <- crossprod(influence[,(q+1):(q+k),drop=FALSE])/n^2+
      sensitivity%*%Vm%*%t(sensitivity)
    close_to(f$vcov,expected,1e-7)
  }
  rows <- maic_wald_rows(z)
  stopifnot(nrow(rows)==18L,all(subset(rows,method %in% methods)$status=="available"),
    all(subset(rows,method %in% methods)$fitting_weights==if(mode=="unit")"ones" else "MAIC"))
  if (mode=="unit") for(method in methods) {
    f <- z$fits[[method]]; d0 <- d1 <- dat; d0$x <- 0; d1$x <- 1
    X0 <- model.matrix(delete.response(terms(f$fit)),d0)
    X1 <- model.matrix(delete.response(terms(f$fit)),d1)
    p0 <- predict(f$fit,d0,type="response"); p1 <- predict(f$fit,d1,type="response")
    v <- z$weights$weight; m0 <- weighted.mean(p0,v); m1 <- weighted.mean(p1,v)
    gradient <- colSums(v*(X1*p1/m1-X0*p0/m0))/sum(v)
    expected <- sqrt(drop(t(gradient)%*%f$vcov%*%gradient))
    rr <- rows[rows$method==method & rows$ci_type=="Wald",]
    close_to(rr$estimate,log(m1/m0)); close_to(rr$se,expected)
  }
  batch <- suppressWarnings(simulate_anchored(design,500,seeds=811,score_spec=spec,
    weight_covariates=paste0("z",1:6),model_weighting=mode))
  stopifnot(all(methods %in% batch$reference_truth$method),
    all(methods %in% summarize_anchored_simulation(batch)$metrics$method),
    length(batch$diagnostics[[1]]$glm_fc)==6L)
  cat("PASS:",mode,"coefficient fits, weights, joint covariance, marginal Wald and simulation truth\n")
}

## A failed new method retains the schema and has a separate last-iterate archive.
failed <- suppressWarnings(run.anchored(design,100,seed=811,score_spec=spec,
  weight_covariates=paste0("z",1:6),fit_control=list(fc_max_step=1L,fc_thres=1e-14)))
stopifnot(all(methods %in% failed$estimates$method),
  any(methods %in% names(failed$nonconverged_fits)),
  all(is.na(subset(failed$estimates,method %in% methods & status!="available")$estimate)))
cat("PASS: failed FC methods retain named rows and archived last iterates\n")
