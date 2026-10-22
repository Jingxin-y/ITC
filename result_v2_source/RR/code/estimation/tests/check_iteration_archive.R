## Explicit MAIC mode preserves the weighted regression fixtures; unit mode has its own checks.
## Run from anchored_v1: Rscript estimation/tests/check_iteration_archive.R
source("load_anchored.R",encoding="UTF-8"); source_anchored()
close_to <- function(a,b,tol=1e-8) stopifnot(isTRUE(all.equal(a,b,tolerance=tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
limits <- anchored_fit_control()
stopifnot(identical(limits,list(brm_max_step=1000L,brm_optim_maxit=200L,
  fc_max_step=300L,glm_maxit=200L,brm_thres=1e-6,fc_thres=1e-6,brm_optim_reltol=1e-6)))
for (x in list(list(20),list(unknown=20),list(glm_maxit=0),list(glm_maxit=1.5),
               list(fc_max_step=NA_real_),list(brm_max_step=NULL))) bad(anchored_fit_control(x))
design <- joint_rr_gop_design(event="common",ess_ratio=.5)
scores <- make_brm_score_spec()
good <- suppressWarnings(run.anchored(model_weighting="maic",design,100,100,seed=2,score_spec=scores))
stopifnot(isTRUE(good$fits$LB$converged),good$fits$LB$control$maxit==200L,
  good$fits$brm$optimization$optim_maxit==200L,
  good$fits$brm$optimization$optim_reltol==1e-6,
  good$fits$brm$optimization$max_step==1000L,
  identical(good$fits[["brm-FC"]]$fit_control,limits),
  nrow(subset(good$nonconverged_estimates,!method %in% anchored_glm_fc_method_names()))==0L,
  all(names(good$nonconverged_fits) %in% anchored_glm_fc_method_names()))
cat("PASS: increased limits reach each solver; the previously slow LB fit now converges\n")

short <- anchored_fit_control(list(brm_max_step=1L,brm_optim_maxit=1L,fc_max_step=1L,glm_maxit=1L))
checkpoint <- tempfile(fileext=".rds")
batch <- suppressWarnings(simulate_anchored(model_weighting="maic",design,100,100,seeds=2,score_spec=scores,
  fit_control=short,checkpoint_file=checkpoint))
raw <- batch$nonconverged_estimates
stopifnot(nrow(raw)>0,all(raw$r==2),all(raw$status=="fit_not_converged"),
  any(is.finite(raw$estimate)),all(is.na(subset(batch$estimates,status=="fit_not_converged")$estimate)),
  identical(readRDS(checkpoint)$nonconverged_estimates,raw))
fits <- readRDS(checkpoint)$nonconverged_fits[[1]]$fits
for (method in c("brm","brm-FC")) {
  f <- fits[[method]]
  stopifnot(isFALSE(f$converged),f$last_optimizer$step==1L,
    all(is.finite(f$point.est)),all(is.na(f$se.est)),identical(f$fit_control,short))
  close_to(unname(c(f$optimization_input$alpha.start,f$optimization_input$beta.start)),
           unname(f$point.est))
  retry <- f$optimization_input
  if (method=="brm") stopifnot(retry$optim.reltol==1e-6)
  retry$max.step <- if (method=="brm") 200L else 300L
  if (method=="brm") retry$optim.maxit <- 200L
  solver <- if (method=="brm") max.likelihood else max.likelihood.firth.rr
  resumed <- suppressWarnings(do.call(solver,retry))
  stopifnot(isTRUE(resumed$convergence))
}
stopifnot(fits$LB$control$maxit==1L,fits$LP$control$maxit==1L,
          fits$RLP$fit$control$maxit==1L)
metrics <- summarize_anchored_simulation(batch)$metrics
stopifnot(all(subset(metrics,method %in% c("brm","brm-FC","LB"))$n_estimate==0))
cat("PASS: last BRM/FC/GLM values survive checkpointing separately; saved BRM/FC inputs can resume to convergence\n")

## This checks status/coverage plumbing using explicit diagnostic fixtures, not BC accuracy.
bc <- list(low=-2,up=2,p=.25,diagnostics=list(
  lower=list(list(convergence=TRUE)),
  upper=list(list(convergence=FALSE,reason="confirmed rejected point not found in search range")),
  observed.null=list(list(converged=TRUE)),p.value=list(list(diagnostics=list()))))
stopifnot(bc_parameter_status(bc,1)$interval_status=="BC_boundary_truncated")
d <- subset(good$estimates,method=="brm" & term=="(Intercept)" & component=="alpha")
d$method <- "brm-BC"; d$interval <- "BC"; d$lower <- -2; d$upper <- 2
d$interval_status <- bc_parameter_status(bc,1)$interval_status
fixture <- structure(list(estimates=cbind(r=2,d),seeds=2),class=c("anchored_simulation","list"))
truth <- data.frame(method="brm-BC",component="alpha",term="(Intercept)",truth=0)
m <- summarize_anchored_simulation(fixture,truth)$metrics
stopifnot(m$n_interval==1,m$n_boundary_truncated==1,m$n_coverage==1,m$coverage==1)
bc$diagnostics$lower[[1]] <- list(convergence=FALSE,reason="maximum bisection steps reached before width threshold")
stopifnot(bc_parameter_status(bc,1)$interval_status=="BC_search_failed")
bc$diagnostics$lower[[1]] <- list(convergence=TRUE)
bc$diagnostics$upper[[1]]$diagnostics <- list(optimizer.failures=1)
stopifnot(bc_parameter_status(bc,1)$interval_status=="BC_fit_failed")
bc$diagnostics$upper[[1]]$diagnostics <- list()
bc$low <- 3
stopifnot(bc_parameter_status(bc,1)$interval_status=="BC_invalid_interval")
stopifnot(formals(exact)$optim.maxit==300,formals(exact_safe)$max.step==100)
cat("PASS: pure boundary truncation counts in coverage; search/fit failure or invalid endpoints remain excluded\n")

## Actual CLI persistence, including the archive, occurs before summary.
outdir <- tempfile("anchored_archive_"); dir.create(outdir)
config <- file.path(outdir,"config.R")
writeLines(c(
  'simulation_args <- list(design=joint_rr_gop_design(event="common",ess_ratio=.5),',
  ' n_AC=100,n_BC=100,seeds=2,score_spec=make_brm_score_spec(),run_exact=FALSE,model_weighting="maic",',
  ' fit_control=list(brm_max_step=1L,brm_optim_maxit=1L,fc_max_step=1L,glm_maxit=1L))'),config)
log <- system2(file.path(R.home("bin"),"Rscript.exe"),
  c("--vanilla",shQuote("estimation/run_simulation_anchored.R"),shQuote(config),shQuote(outdir)),
  stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(log,"status")))
disk <- readRDS(file.path(outdir,"anchored_joint_nonconverged.rds"))
csv <- read.csv(file.path(outdir,"anchored_joint_nonconverged.csv"))
stopifnot(nrow(csv)==nrow(raw),length(disk$fits)>0,identical(disk$fit_control,short),
          identical(disk$estimates,raw),isFALSE(disk$fits[[1]]$fits[["brm-FC"]]$converged))
cat("PASS: CLI writes separate nonconverged CSV/RDS with last fits and settings\n")
