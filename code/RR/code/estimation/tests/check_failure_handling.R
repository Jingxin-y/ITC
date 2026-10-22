## Explicit MAIC mode preserves the weighted regression fixtures; unit mode has its own checks.
## Run from anchored_v1: Rscript estimation/tests/check_failure_handling.R
source("load_anchored.R", encoding="UTF-8")
source_anchored()
close_to <- function(a,b,tol=1e-8) stopifnot(isTRUE(all.equal(a,b,tolerance=tol)))
capture <- function(expr) tryCatch(force(expr), error=function(e)e)
design <- joint_rr_gop_design(event="common",ess_ratio=.5)
scores <- make_brm_score_spec()
result <- suppressWarnings(run.anchored(model_weighting="maic",score_standardization="fixed",design,100,100,seed=2,score_spec=scores,
  weight_covariates=paste0("z",1:6),return_data=TRUE,fit_control=list(glm_maxit=25L)))
lb <- subset(result$estimates,method=="LB")
stopifnot(!isTRUE(result$fits$LB$converged), all(lb$status=="fit_not_converged"),
          all(is.na(lb$estimate)), all(is.na(lb$lower)))
bad_brm <- result$fits$brm
bad_brm$converged <- FALSE
b <- comparison_parameter_rows("brm",bad_brm,"brm")
stopifnot(all(b$status=="fit_not_converged"),all(is.na(b$estimate)))
cat("PASS: real LB nonconvergence and BRM convergence flag are excluded\n")

err <- capture(get.weight(cbind(z1=rep(0,20)),c(z1=.7)))
stopifnot(inherits(err,"maic_calibration_error"),all(is.na(err$result$weight)),
          err$result$diagnostics$reason=="target_outside_support")
## Both targets individually lie inside the support; the joint target is impossible.
z <- cbind(z1=rep(0:1,10),z2=rep(0:1,10))
err <- capture(get.weight(z,c(z1=.2,z2=.8)))
stopifnot(inherits(err,"maic_calibration_error"),all(is.na(err$result$weight)))
good <- get.weight(cbind(z1=rep(0:1,10)),c(z1=.7))
close_to(weighted.mean(rep(0:1,10),good$weight),.7,tol=1e-6)
stopifnot(good$status=="available",good$diagnostics$max_scaled_error<=1e-6)
bad_ad <- result$data$data.AD
bad_ad$anchored.mean.x["z1"] <- 2
failed_weights <- suppressWarnings(get.compare.anchored(model_weighting="maic",score_standardization="fixed",result$data$data.IPD,bad_ad,
  score_spec=scores,weight_covariates=paste0("z",1:6),run_exact=FALSE))
stopifnot(failed_weights$weights$status=="weighting_failed",
  all(subset(failed_weights$estimates,method=="brm")$status=="weighting_failed"),
  all(subset(failed_weights$estimates,!method %in% c("brmad","brm_b-BC"))$status=="weighting_failed"),
  is.null(failed_weights$fits$brm))
cat("PASS: infeasible MAIC targets rejected, valid balance checked, configured BRM failure reported\n")

pa <- ncol(result$fits$brm$va)
bc <- list(low=rep(-2,pa),up=rep(2,pa),p=rep(.25,pa),diagnostics=list(
  upper=rep(list(list(convergence=TRUE)),pa), lower=rep(list(list(convergence=TRUE)),pa),
  observed.null=rep(list(list(converged=TRUE)),pa),p.value=rep(list(list(diagnostics=list())),pa)))
bc$diagnostics$upper[[1]] <- list(convergence=FALSE,
  reason="confirmed rejected point not found in search range")
bc$diagnostics$lower[[2]] <- list(convergence=FALSE,reason="maximum bisection steps reached before width threshold")
rows <- comparison_parameter_rows("brm-BC",result$fits$brm,"brm",bc,use_bc=TRUE)
alpha <- subset(rows,component=="alpha")
stopifnot(all(alpha$status=="available"),all(is.finite(alpha$estimate)),
  identical(alpha$interval_status,c("BC_boundary_truncated","BC_search_failed")),
  alpha$lower[1]==-2,alpha$upper[1]==2,is.na(alpha$lower[2]),is.na(alpha$upper[2]),
  all(alpha$p_status=="available"))
bc$diagnostics <- NULL
unverified <- comparison_parameter_rows("brm-BC",result$fits$brm,"brm",bc,use_bc=TRUE)
stopifnot(all(subset(unverified,component=="alpha")$interval_status=="BC_unverified"))
cat("PASS: BC boundary truncation retains endpoints; real search failures and missing diagnostics are excluded\n")

## Deliberately finite failed rows verify that summary trusts statuses, not just finiteness.
lb$estimate <- .2; lb$se <- .1; lb$lower <- -.1; lb$upper <- .5; lb$p.value <- .5
alpha$lower <- -2; alpha$upper <- 2
batch <- structure(list(estimates=cbind(r=1,rbind(lb,alpha)),seeds=1),class=c("anchored_simulation","list"))
truth <- data.frame(method=c("LB","brm-BC"),component=c("treatment","alpha"),
  term=c("x","h_E"),truth=c(.2,.3))
s <- summarize_anchored_simulation(batch,truth)
stopifnot(all(subset(s$metrics,method=="LB")$n_estimate==0),
  all(subset(s$metrics,method=="LB")$estimate_failure_rate==1),
  subset(s$metrics,method=="brm-BC" & term=="(Intercept)")$n_interval==1,
  subset(s$metrics,method=="brm-BC" & term=="(Intercept)")$n_boundary_truncated==1,
  subset(s$metrics,method=="brm-BC" & term=="h_E")$n_interval==0,
  all(is.na(s$metrics$coverage)))
cat("PASS: failed estimates/searches excluded; boundary-truncated intervals counted separately\n")

## Exceptions returning NULL must also preserve the complete parameter schema.
fits <- result$fits
fits["brm-FC"] <- list(NULL)
schemas <- comparison_schema(apply_brm_scores(result$data$data.IPD,scores)$data,
  c("h_E","h_P"),"h_E",list(va=result$fits$brm$va,vb=result$fits$brm$vb))
missing <- finish_anchored_comparison(fits,result$bc,result$weights,result$variables,
                                     schemas,FALSE,TRUE)
missing_batch <- structure(list(estimates=cbind(r=2,missing$estimates),seeds=2,
  parameter_schema=missing$parameter_schema),class=c("anchored_simulation","list"))
fc_truth <- data.frame(method="brm-FC",component="alpha",term="h_E",truth=.3)
fc_summary <- summarize_anchored_simulation(missing_batch,fc_truth)$metrics
stopifnot(all(subset(fc_summary,method=="brm-FC")$n_estimate==0),
          subset(fc_summary,method=="brm-FC" & term=="h_E")$truth==.3)
cat("PASS: NULL fit retains named parameters and accepts correct truth for an entirely failed method\n")

edge <- result$data$data.IPD$data
edge$y[edge$x==0] <- 0
set.seed(1)
fallback <- get.brm.ad(make_na_est(),edge,rep(1,nrow(edge)))
fallback_rows <- comparison_parameter_rows("brmad",fallback,"brm")
stopifnot(fallback_rows$status=="available",is.finite(fallback_rows$estimate),
          fallback_rows$interval=="Beta posterior",fallback_rows$interval_status=="available")
cat("PASS: scalar Bayes fallback retains its valid point/interval when ordinary BRM has failed\n")

## A failure on replicate 2 must leave replicate 1 readable on disk.
real_run <- run.anchored
checkpoint <- tempfile(fileext=".rds")
run.anchored <- function(...) {
  a <- list(...)
  if (a$seed==3) stop("injected replicate failure")
  do.call(real_run,a)
}
err <- tryCatch(suppressWarnings(capture(simulate_anchored(model_weighting="maic",score_standardization="fixed",design,100,100,seeds=2:3,
  score_spec=scores,weight_covariates=paste0("z",1:6),checkpoint_file=checkpoint))),
  finally={run.anchored <- real_run})
saved <- readRDS(checkpoint)
stopifnot(inherits(err,"error"),identical(saved$seeds,2L),!saved$complete,
          saved$interrupted_at_seed==3,nrow(saved$estimates)>0)
cat("PASS: unexpected replicate failure preserves completed work and records interrupted seed\n")

## A CLI summary error occurs after raw RDS/CSV persistence.
outdir <- tempfile("anchored_cli_"); dir.create(outdir)
config <- file.path(outdir,"config.R")
writeLines(c(
  'simulation_args <- list(design=joint_rr_gop_design(event="common",ess_ratio=.5),',
  ' n_AC=100,n_BC=100,seeds=2,score_spec=make_brm_score_spec(),run_exact=FALSE)',
  'truth <- data.frame(method="unknown_method",component="alpha",term="h_E",truth=.3)'),config)
log <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript.exe"),
  c("--vanilla",shQuote("estimation/run_simulation_anchored.R"),shQuote(config),shQuote(outdir)),
  stdout=TRUE,stderr=TRUE))
stopifnot(!is.null(attr(log,"status")),file.exists(file.path(outdir,"anchored_joint.rds")),
          file.exists(file.path(outdir,"anchored_joint_estimates.csv")),
          readRDS(file.path(outdir,"anchored_joint.rds"))$complete)
cat("PASS: summary failure leaves completed CLI raw RDS and estimates CSV intact\n")
