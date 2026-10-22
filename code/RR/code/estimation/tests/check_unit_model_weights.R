## Model weights are one; target standardization must still use MAIC weights.
source("load_anchored.R",encoding="UTF-8")
source_anchored()
close_to <- function(a,b,tol=1e-8)
  stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))

design <- webster_rr_design(ess_ratio=.5)
scores <- make_brm_score_spec()
result <- suppressWarnings(run.anchored(design,1000,1200,seed=811,
  score_spec=scores,weight_covariates=paste0("z",1:6),return_data=TRUE))
ones <- rep(1,1000)
stopifnot(result$variables$model_weighting=="unit",
  result$variables$standardization_weighting=="maic",
  all(result$model_weights$weight==1),length(result$model_weights$ga)==0L,
  result$weights$status=="available",sd(result$weights$weight)>.1)
ipd <- apply_brm_scores(result$data$data.IPD,result$variables$score_spec)
dat <- ipd$data
close_to(colSums(as.matrix(dat[paste0("z",1:6)])*result$weights$weight)/sum(result$weights$weight),
         result$data$data.AD$anchored.mean.x,tol=1e-6)
manual_brm <- get.estimate("AD",ipd,ones,score_spec=result$variables$score_spec)
close_to(result$fits$brm$point.est,manual_brm$point.est)
close_to(result$fits$brm$cov,manual_brm$cov)
stopifnot(!anyDuplicated(names(result$fits)), !"brm_an" %in% names(result$fits))
fc <- get.estimate("AD",ipd,ones,method="firth",score_spec=result$variables$score_spec)
close_to(result$fits[["brm-FC"]]$point.est,fc$point.est)
close_to(result$fits[["brm-FC"]]$cov,fc$cov)
stopifnot(result$fits[["brm-FC"]]$variance.method=="fixed_weight_model_sandwich_first_order")
cat("PASS: calibration still runs; BRM and FC use fixed unit weights and fixed-weight covariance\n")

form <- y~x+h_E+h_P+x:h_E
poisson <- glm(form,data=dat,family=poisson(link="log"))
binomial <- glm(form,data=dat,family=binomial(link="log"),
  start=log_link_start(dat,model_covariates=c("h_E","h_P"),effect_modifiers="h_E"),
  control=glm.control(maxit=200))
close_to(coef(result$fits$LP),coef(poisson))
close_to(coef(result$fits$LB),coef(binomial))
for (method in c("LB","LP")) stopifnot(all(result$fits[[method]]$prior.weights==1))
## Independently assemble HC0 covariance, with no calibration nuisance scores.
for (method in c("RLP","RLB")) {
  f <- result$fits[[method]]$fit
  X <- model.matrix(f); mu <- fitted(f); y <- dat$y
  is_lb <- method=="RLB"
  residual <- if (is_lb) (y-mu)/(1-mu) else y-mu
  curvature <- if (is_lb) mu*(1-y)/(1-mu)^2 else mu
  bread <- MASS::ginv(crossprod(X,X*curvature))
  expected <- bread %*% crossprod(X*residual) %*% t(bread)
  close_to(result$fits[[method]]$vcov,expected)
  stopifnot(all(result$fits[[method]]$vcov.target==0),all(f$prior.weights==1))
}
cmh <- get.weighted.RR(dat,result$weights$weight,model_covariates=c("h_E","h_P"))
close_to(result$fits$CMH,cmh)
close_to(result$fits$CMH_sandwich[c("est","se","low","up","p")],
         cmh[c("est","se","low","up","p")])
stopifnot(result$variables$cmh_weighting=="maic",
          result$variables$cmh_sandwich_uncertainty=="fixed")
cat("PASS: GLM fits use ones; CMH uses MAIC counts with fixed-weight variance\n")

## Different calibration choices must not change model fits or coefficient SEs.
other <- suppressWarnings(run.anchored(design,1000,1200,seed=811,
  score_spec=scores,weight_covariates="z4"))
a <- subset(result$estimates,component!="marginal" & !method %in% c("CMH","CMH_sandwich"))
b <- subset(other$estimates,component!="marginal" & !method %in% c("CMH","CMH_sandwich"))
close_to(a,b)
stopifnot(max(abs(result$weights$weight-other$weights$weight))>.1)
close_to(other$fits$CMH,get.weighted.RR(dat,other$weights$weight,
                                      model_covariates=c("h_E","h_P")))
stopifnot(abs(result$fits$CMH$est-other$fits$CMH$est)>1e-3)
for (method in c("brm-BC","brm-FC-BC")) {
  f <- result$fits[[method]]
  expected <- marginal_logrr_target(f$point.est,f$va,f$vb,result$weights$weight)
  unweighted <- marginal_logrr_target(f$point.est,f$va,f$vb,ones)
  actual <- result$marginal_estimates$estimate[result$marginal_estimates$method==method]
  close_to(actual,expected)
  stopifnot(abs(actual-unweighted)>1e-3)
}
cat("PASS: fitting ignores calibration choices while marginal effects use BC target weights\n")

## Check both bootstrap entry points, including their nuisance-model weights.
real_marginal <- exact_marginal_target
calls <- list()
exact_marginal_target <- function(...) {
  a <- list(...); calls[[length(calls)+1L]] <<- a
  p <- marginal_logrr_target(a$pars,a$va,a$vb,a$standardization.weight)
  list(estimate=p,low=p-.5,up=p+.5,p=.5,status="available",interval.available=TRUE,
    p.status="available",diagnostics=list(marker="weight separation"))
}
bc <- tryCatch(suppressWarnings(run.anchored(design,1000,1200,seed=811,
  score_spec=scores,weight_covariates=paste0("z",1:6),run_exact=TRUE)),
  finally={exact_marginal_target <- real_marginal})
stopifnot(length(calls)==2L)
for (a in calls) {
  close_to(a$weight,ones)
  close_to(a$standardization.weight,bc$weights$weight)
}
real_exact <- exact; coefficient_calls <- list()
exact <- function(...) {
  a <- list(...); coefficient_calls[[length(coefficient_calls)+1L]] <<- a
  list(low=rep(-2,a$pa),up=rep(2,a$pa),p=rep(.5,a$pa))
}
tryCatch(suppressWarnings(run.anchored(design,1000,1200,seed=811,
  score_spec=scores,weight_covariates=paste0("z",1:6),run_exact=TRUE,bc_target="coefficient")),
  finally={exact <- real_exact}) -> coefficient_result
stopifnot(length(coefficient_calls)==2L,all(vapply(coefficient_calls,
  function(a) all(a$weight==1),logical(1))))
cat("PASS: marginal BC separates fitting/target weights; coefficient BC fits with ones\n")

## Calibration failure blocks target effects, not independent unweighted fits.
bad_ad <- result$data$data.AD; bad_ad$anchored.mean.x["z1"] <- 2
failed <- suppressWarnings(get.compare.anchored(result$data$data.IPD,bad_ad,
  score_spec=scores,score_standardization="fixed",weight_covariates=paste0("z",1:6),run_exact=FALSE))
stopifnot(failed$weights$status=="weighting_failed",failed$model_weights$status=="available",
  all(subset(failed$estimates,method %in% c("brm","brm-FC","LP"))$status=="available"),
  all(subset(failed$estimates,method %in% c("CMH","CMH_sandwich"))$status=="weighting_failed"),
  is.null(failed$fits$CMH),is.null(failed$fits$CMH_sandwich),
  all(subset(failed$marginal_estimates,method!="brm_b-BC")$status=="weighting_failed"))
cat("PASS: calibration failure retains unweighted fits and marks unavailable target effects\n")

## Preserve explicit access to the former MAIC-weighted behavior.
weighted <- suppressWarnings(run.anchored(design,1000,1200,seed=811,
  score_spec=scores,weight_covariates=paste0("z",1:6),model_weighting="maic"))
close_to(weighted$model_weights$weight,weighted$weights$weight)
close_to(weighted$fits$CMH,result$fits$CMH)
joint_cmh <- get.weighted.RR.sandwich(dat,weighted$weights$weight,
  gamma=weighted$weights$ga,target_mean=result$data$data.AD$anchored.mean.x,
  target_cov=result$data$data.AD$anchored.cov.x,m_AD=result$data$data.AD$anchored.n,
  account_for_gamma=TRUE,account_for_target_mean=TRUE,
  weight_covariates=weighted$weights$covariates,model_covariates=c("h_E","h_P"))
close_to(weighted$fits$CMH_sandwich,joint_cmh)
stopifnot(weighted$variables$cmh_sandwich_uncertainty=="joint")
manual <- glm(form,data=dat,family=poisson(link="log"),weights=weighted$weights$weight)
## GLM starting values differ; both solutions should agree to solver accuracy.
close_to(coef(weighted$fits$LP),coef(manual),tol=1e-6)
stopifnot(max(abs(coef(weighted$fits$LP)-coef(result$fits$LP)))>1e-4)
batch <- suppressWarnings(simulate_anchored(design,1000,1200,seeds=811,
  score_spec=scores,weight_covariates=paste0("z",1:6)))
stopifnot(batch$model_weighting=="unit",batch$diagnostics[[1]]$model_ess==1000,
          batch$diagnostics[[1]]$ess<1000)
close_to(subset(batch$estimates,select=-r),result$estimates)
cat("PASS: explicit weighted compatibility and simulation diagnostics distinguish both ESS values\n")
