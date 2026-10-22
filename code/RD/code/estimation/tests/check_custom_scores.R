## Coefficient-BC compatibility; marginal integration is checked in check_marginal_pipeline.R.
## Run from anchored_v1: Rscript estimation/tests/check_custom_scores.R
source("load_anchored.R", encoding="UTF-8")
source_anchored()
close_to <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(a,b,tolerance=tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
zcols <- paste0("z",1:6)
z <- expand.grid(setNames(rep(list(0:1),6),zcols))
default <- make_brm_score_spec()
scored <- apply_brm_scores(z,default)
original <- rr_gop_scores(z)
for (name in c("S_E","S_P","h_E","h_P")) close_to(scored[[name]],original[[name]])
stopifnot(length(unique(scored$h_E))==4L,length(unique(scored$h_P))==4L)
close_to(unname(default$center),c(.9,.9))
close_to(unname(default$scale),rep(sqrt(.63),2))

spec <- make_brm_score_spec(em_covariates=c("z1","z4"),
                           non_em_covariates=c("z2","z3","z5","z6"))
custom <- apply_brm_scores(z,spec)
close_to(custom$S_E,1-z$z1+z$z4)
close_to(custom$h_E,(1-z$z1+z$z4-.6)/sqrt(.42))
close_to(custom$h_P,(z$z2+1-z$z3+1-z$z5+z$z6-1.2)/sqrt(.84))
stopifnot(length(unique(custom$h_E))==3L,length(unique(custom$h_P))==5L)
close_to(apply_brm_scores(z[1:8,],spec)$h_E,custom$h_E[1:8])
weighted_spec <- make_brm_score_spec(c("z1","z4"),character(0),
  reverse_covariates=character(0),coefficients=c(z1=2,z4=-1))
close_to(apply_brm_scores(z,weighted_spec)$h_E,
         (2*z$z1-z$z4-1.1)/sqrt(1.05))
cat("PASS: default and custom scores, directions, coefficients and fixed standardization\n")

## Balanced observations, non-boundary events in every treatment/covariate cell.
dat <- merge(z,expand.grid(x=0:1,id=1:12),by=NULL)
dat$y <- as.integer(dat$id <= 2+dat$x+dat$z1+dat$z4)
for (name in c("S_E","S_P","h_E","h_P")) dat[[name]] <- rr_gop_scores(dat)[[name]]
ipd <- list(data=dat,truth=list(unchanged=TRUE))
saved <- ipd
set.seed(42)
ad_individual <- list(data=dat[sample(nrow(dat),700,replace=TRUE,
  prob=1+dat$z1*dat$z2+dat$z4),])
ad <- get.aggre(ad_individual)
saved_ad <- ad
score_cols <- c("h_E","h_P","S_E","S_P")
ad_auto <- brm_score_summary(ad,spec,score_cols)
ad_manual <- get.aggre(apply_brm_scores(ad_individual,spec),c(zcols,score_cols))
close_to(ad_auto,ad_manual)
stopifnot(identical(ad,saved_ad))
bad(brm_score_summary(get.aggre(ad_individual,"h_E"),spec,"h_E"))
cat("PASS: AD score means/covariance match individual-data aggregation, including correlations\n")

wf <- get.weight(covariate_matrix(ipd$data,zcols),ad$anchored.mean.x)
args <- list(group="IPD",data=ipd,weight=wf$weight,gamma.ipd=wf$ga,
  target_mean=ad$anchored.mean.x,target_cov=ad$anchored.cov.x,m_AD=ad$anchored.n,
  weight_covariates=zcols,score_spec=spec)
fit <- do.call(get.estimate,args)
manual <- ipd
manual$data$S_E <- 1+dat$z4-dat$z1
manual$data$S_P <- dat$z2+1-dat$z3+1-dat$z5+dat$z6
manual$data$h_E <- (manual$data$S_E-.6)/sqrt(.42)
manual$data$h_P <- (manual$data$S_P-1.2)/sqrt(.84)
manual_args <- args
manual_args$data <- manual
manual_args$score_spec <- NULL
manual_args$effect_modifiers <- "h_E"
manual_args$non_effect_modifiers <- "h_P"
manual_args$model_covariates <- c("h_E","h_P")
manual_fit <- do.call(get.estimate,manual_args)
close_to(fit$point.est,manual_fit$point.est)
close_to(fit$cov,manual_fit$cov)
stopifnot(ncol(fit$va)==2L,ncol(fit$vb)==3L,
          identical(fit$score_spec,spec),identical(ipd,saved))
prepared_args <- args
prepared_args$data <- apply_brm_scores(ipd,spec)
prepared_args$score_spec <- NULL
close_to(do.call(get.estimate,prepared_args)$point.est,fit$point.est)
cat("PASS: BRM estimates and joint covariance match manually constructed custom scores\n")

## Record BC forwarding only; no claim about interval accuracy/coverage.
real_exact <- exact
recorded <- list()
exact <- function(...) {
  a <- list(...)
  recorded[[length(recorded)+1L]] <<- a
  list(low=rep(-1,a$pa),up=rep(1,a$pa),p=rep(.5,a$pa))
}
bc <- tryCatch(get.brm.bc(bc_target="coefficient",ipd,ad,wf$weight,wf$ga,fit,NULL),
               finally={exact <- real_exact})
stopifnot(length(recorded)==2L,identical(bc$firth$score_spec,spec),
          identical(bc$firth$va,fit$va),identical(bc$firth$vb,fit$vb),
          identical(recorded[[1]]$va,fit$va),identical(recorded[[2]]$vb,bc$firth$vb))
firth_args <- manual_args
firth_args$method <- "firth"
close_to(bc$firth$point.est,do.call(get.estimate,firth_args)$point.est)
cat("PASS: Firth refit restores score specification and BC uses the fitted matrices\n")

comparison <- suppressWarnings(get.compare.anchored(bc_target="coefficient",score_standardization="fixed",ipd,ad,run_exact=FALSE,
  score_spec=spec,weight_covariates="h_E"))
manual_comparison <- suppressWarnings(get.compare.anchored(bc_target="coefficient",score_standardization="fixed",manual,ad_manual,
  run_exact=FALSE,weight_covariates="h_E",model_covariates=c("h_E","h_P"),
  effect_modifiers="h_E",non_effect_modifiers="h_P"))
close_to(comparison$weights[setdiff(names(comparison$weights),"diagnostics")],
         manual_comparison$weights[setdiff(names(manual_comparison$weights),"diagnostics")])
## Near-zero balance residuals need an absolute tolerance, not relative error.
stopifnot(isTRUE(all.equal(comparison$weights$diagnostics,manual_comparison$weights$diagnostics,
                          tolerance=1e-9,scale=1)))
close_to(comparison$estimates,manual_comparison$estimates,tol=1e-7)
available_rows <- !comparison$estimates$method %in% c("brmad","brm_b-BC")
stopifnot(all(comparison$estimates$status[available_rows]=="available"),
          all(is.finite(comparison$estimates$estimate[available_rows])))
stopifnot(identical(comparison$variables$score_spec,spec),
          identical(ipd,saved),identical(ad,saved_ad))
score_bc <- get.brm.bc(bc_target="coefficient",ipd,ad,
  comparison$model_weights$weight,comparison$model_weights$ga,
  comparison$fits$brm,NULL,run_exact=FALSE,
  weight_covariates=comparison$model_weights$covariates,
  model_weighting=comparison$variables$model_weighting,
  standardization.weight=comparison$weights$weight)
close_to(score_bc$firth$point.est,comparison$fits[["brm-FC-BC"]]$point.est)
cat("PASS: all comparison methods and score weighting match manual configuration\n")

none <- make_brm_score_spec(character(0),character(0))
empty_args <- args
empty_args$score_spec <- none
empty <- do.call(get.estimate,empty_args)
stopifnot(ncol(empty$va)==1L,ncol(empty$vb)==1L)
args$effect_modifiers <- character(0)
disabled <- do.call(get.estimate,args)
stopifnot(ncol(disabled$va)==1L,ncol(disabled$vb)==1L)
bad(make_brm_score_spec("z1","z1"))
bad(make_brm_score_spec("z1",character(0),reverse_covariates="z2"))
bad(make_brm_score_spec("unknown",character(0)))
bad(make_brm_score_spec("y",character(0)))
bad(make_brm_score_spec("z1",character(0),reference_prob=c(z1=0)))
bad(apply_brm_scores(z[,-1],default))
not_binary <- z
not_binary$z1[1] <- .2
bad(apply_brm_scores(not_binary,default))
cat("PASS: empty scores, explicit no-EM override and invalid configuration rejection\n")

design <- joint_rr_gop_design(.4,.2,.3,.5,-3,.2,.4,ess_ratio=.5)
set.seed(12)
sim <- generate_anchored_data(100,100,design)
original_sim <- sim
prepared_sim <- apply_brm_scores(sim$data.IPD,spec)
stopifnot(identical(sim,original_sim),
  identical(prepared_sim$truth,sim$data.IPD$truth),
  identical(prepared_sim$design,sim$data.IPD$design),
  identical(prepared_sim$data$y,sim$data.IPD$data$y),
  identical(prepared_sim$data$x,sim$data.IPD$data$x),
  identical(prepared_sim$data[zcols],sim$data.IPD$data[zcols]))
cat("PASS: DGM data, risks and truth remain unchanged\n")
