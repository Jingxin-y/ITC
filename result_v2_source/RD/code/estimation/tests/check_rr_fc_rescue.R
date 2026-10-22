## Run from anchored_v1:
## Rscript estimation/tests/check_rr_fc_rescue.R . failed_fc_details.rds output_dir
## Optional arguments 4/5: pre-merge RR source; historical checkpoint path list.
args <- commandArgs(trailingOnly=TRUE)
code <- normalizePath(if(length(args))args[1L]else".",winslash="/",mustWork=TRUE)
input <- if(length(args)>=2L)args[2L]else file.path(code,
  "result_cluster/result_v2/plots/all_20261008_103434/failed_fc_details.rds")
input <- normalizePath(input,winslash="/",mustWork=TRUE)
out <- if(length(args)>=3L)args[3L]else file.path(tempdir(),"rr_fc_rescue_validation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,winslash="/",mustWork=TRUE)
if(.Platform$OS.type=="windows") {
  invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
}
stopifnot(requireNamespace("nleqslv",quietly=TRUE))
source(file.path(code,"load_anchored.R"),encoding="UTF-8");source_anchored(code)
merged.rr <- max.likelihood.firth.rr
helper <- .firth_rr_rescue
if(length(args)>=4L) {
  old.env <- new.env(parent=.GlobalEnv)
  sys.source(args[4L],envir=old.env)
  baseline.rr <- old.env$max.likelihood.firth.rr
  stopifnot(identical(body(firth_rr_derivatives),body(old.env$firth_rr_derivatives)))
} else {
  # The unchanged shared core is precisely the old RR public entry point.
  baseline.rr <- function(param,y,x,va,vb,alpha.start,beta.start,weight,max.step,thres,pa,pb)
    max.likelihood.firth(y,x,va,vb,alpha.start,beta.start,weight,max.step,thres,pa,pb,firth_rr_derivatives)
}
with.solver <- function(solver,expr) {
  old <- max.likelihood.firth.rr
  assign("max.likelihood.firth.rr",solver,envir=.GlobalEnv)
  on.exit(assign("max.likelihood.firth.rr",old,envir=.GlobalEnv))
  force(expr)
}
wrapped <- function(a) do.call(MLEst,c(a,list(group="AD",gamma.weight=0,
  target_mean=0,method="firth")))

records <- readRDS(input); original.hash <- tools::md5sum(input)
stopifnot(length(records)==8L)
recovered <- failures <- list()
for(i in seq_along(records)) {
  z <- records[[i]];a <- z$fit$optimization_input
  old <- do.call(baseline.rr,a)
  stopifnot(!old$convergence)
  rng <- get0(".Random.seed",envir=.GlobalEnv)
  fit <- wrapped(a)
  stopifnot(fit$converged,fit$fc.diagnostics$rescue$accepted,
    identical(rng,get0(".Random.seed",envir=.GlobalEnv)),
    fit$fc.diagnostics$iterations<=a$max.step,
    all(is.finite(c(fit$point.est,fit$se.est,fit$cov))))
  terms <- firth_rr_derivatives(fit$point.est[seq_len(a$pa)],
    fit$point.est[a$pa+seq_len(a$pb)],a$x,a$va,a$vb)
  w <- a$weight/max(a$weight);w <- w/mean(w)
  check <- firth_score_adjustment(a$y,a$va,a$vb,w,terms)
  original.score <- max(abs(check$adjusted.score/z$fit$fc.diagnostics$score_scale))
  entry.score <- max(abs(check$adjusted.score/fit$fc.diagnostics$score_scale))
  prob <- getProbRR(drop(a$va%*%fit$point.est[seq_len(a$pa)]),
    drop(a$vb%*%fit$point.est[a$pa+seq_len(a$pb)]))
  stopifnot(terms$interior,original.score<=a$thres,entry.score<=a$thres,
    min(prob)>1e-8,max(prob)<1-1e-8,
    min(eigen(check$information,symmetric=TRUE,only.values=TRUE)$values)>0,
    min(eigen(fit$cov,symmetric=TRUE,only.values=TRUE)$values)>0,
    max(abs(fit$cov-check$covariance))<1e-10)
  # Independently resume the actual archived cluster failure, without rerunning
  # the primary path. Both routes must reach the same root on the same binary y.
  resumed <- helper(z$fit$last_optimizer,a$y,a$x,a$va,a$vb,a$weight,
    a$max.step,a$thres,a$pa,a$pb)
  stopifnot(resumed$convergence,max(abs(fit$point.est-resumed$par))<1e-6)
  failures[[i]] <- data.frame(run=z$run,seed=z$seed,old_converged=old$convergence,
    new_converged=fit$converged,original_scaled_score=original.score,
    entry_scaled_score=entry.score,total_iterations=fit$fc.diagnostics$iterations,
    rescue_iterations=fit$fc.diagnostics$rescue$iterations,
    min_probability=min(prob),max_probability=max(prob))
  recovered[[i]] <- list(run=z$run,seed=z$seed,fit=fit,original_input=a,
    original_rds_md5=unname(original.hash))
  cat("PASS failure",i,"seed",z$seed,"score",original.score,"\n");flush.console()
}
write.csv(do.call(rbind,failures),file.path(out,"eight_failures_validation.csv"),row.names=FALSE)
saveRDS(recovered,file.path(out,"recovered_eight_fits.rds"),version=2)

## Previously successful real checkpoints + common/rare weighted RR samples.
cases <- list()
if(length(args)>=5L) {
  for(path in readLines(args[5L],warn=FALSE,encoding="UTF-8")) {
    z <- readRDS(path);r <- z$result;f <- r$fits[["brm-FC"]]
    stopifnot(isTRUE(f$converged))
    a <- list(param="RR",y=r$data$data.IPD$data$y,x=r$data$data.IPD$data$x,
      va=f$va,vb=f$vb,alpha.start=rep(0,ncol(f$va)),beta.start=rep(0,ncol(f$vb)),
      weight=r$model_weights$weight,max.step=f$fit_control$fc_max_step,
      thres=f$fit_control$fc_thres,pa=ncol(f$va),pb=ncol(f$vb))
    cases[[length(cases)+1L]] <- list(label=path,a=a,archived=f)
  }
}
set.seed(20261019)
for(event in c("common","rare")) for(pa in 1:3) {
  dat <- expand.grid(k=seq_len(if(event=="common")40L else 80L),x=0:1,z1=0:1,z2=0:1)
  va <- cbind(1,dat$z1,dat$z2)[,seq_len(pa),drop=FALSE]
  vb <- cbind(1,dat$z1,dat$z2)
  prob <- getProbRR(drop(va%*%c(.2,.25,-.15)[seq_len(pa)]),
    drop(vb%*%c(if(event=="common")-2 else -5,.3,-.2)))
  y <- rbinom(nrow(dat),1,prob[cbind(seq_len(nrow(dat)),dat$x+1L)])
  for(mode in c("unit","fixed_weight")) {
    a <- list(param="RR",y=y,x=dat$x,va=va,vb=vb,alpha.start=rep(0,pa),
      beta.start=rep(0,3),weight=if(mode=="unit")rep(1,nrow(dat))else exp(.3*dat$z1-.2*dat$z2),
      max.step=300L,thres=1e-6,pa=pa,pb=3L)
    cases[[length(cases)+1L]] <- list(label=paste(event,pa,mode,sep="/"),a=a)
  }
}
unchanged <- list()
# A successful call must never even enter the rescue helper.
assign(".firth_rr_rescue",function(...)stop("Rescue must not run for a converged primary fit"),envir=.GlobalEnv)
for(i in seq_along(cases)) {
  c <- cases[[i]];a <- c$a
  before <- do.call(baseline.rr,a)
  stopifnot(before$convergence)
  before.wrapped <- with.solver(baseline.rr,wrapped(a))
  rng <- .Random.seed
  after <- do.call(merged.rr,a)
  after.wrapped <- wrapped(a)
  archive.delta <- if(is.null(c$archived)) NA_real_ else
    max(abs(before.wrapped$point.est-c$archived$point.est))
  stopifnot(identical(before,after),identical(before.wrapped,after.wrapped),
    identical(rng,.Random.seed),is.null(after$diagnostics$rescue))
  if(!is.null(c$archived)) stopifnot(archive.delta<1e-10)
  unchanged[[i]] <- data.frame(case=c$label,historical=!is.null(c$archived),
    raw_result_identical=identical(before,after),full_MLEst_identical=identical(before.wrapped,after.wrapped),
    rng_identical=identical(rng,.Random.seed),rescue_not_called=TRUE,
    archived_coefficient_max_abs_difference=archive.delta)
}
assign(".firth_rr_rescue",helper,envir=.GlobalEnv)
write.csv(do.call(rbind,unchanged),file.path(out,"successful_results_unchanged.csv"),row.names=FALSE)
cat("PASS",length(cases),"previously converged fits: complete objects identical; rescue never called\n")
flush.console()

## Strict budget, unavailable optional dependency, and weight invariance.
a <- cases[[length(cases)]]$a;a$max.step <- 1L
short <- do.call(merged.rr,a);short.old <- do.call(baseline.rr,a)
stopifnot(!short$convergence,short$step==1L,!short$diagnostics$rescue$attempted,
  identical(short$par,short.old$par))
unavailable <- helper
local.env <- new.env(parent=environment(unavailable))
local.env$requireNamespace <- function(...)FALSE
environment(unavailable) <- local.env
z <- records[[1L]];a <- z$fit$optimization_input
missing <- unavailable(z$fit$last_optimizer,a$y,a$x,a$va,a$vb,a$weight,a$max.step,a$thres,a$pa,a$pb)
stopifnot(!missing$convergence,!missing$diagnostics$rescue$attempted,
  grepl("unavailable",missing$diagnostics$rescue$reason),
  identical(missing$par,z$fit$last_optimizer$par))
scaled <- helper(z$fit$last_optimizer,a$y,a$x,a$va,a$vb,19*a$weight,a$max.step,a$thres,a$pa,a$pb)
normal <- helper(z$fit$last_optimizer,a$y,a$x,a$va,a$vb,a$weight,a$max.step,a$thres,a$pa,a$pb)
stopifnot(scaled$convergence,identical(scaled$par,normal$par),
  identical(tools::md5sum(input),original.hash))
write.csv(data.frame(check=c("strict one-step budget","missing nleqslv retains original failure",
  "global weight scaling invariant","original RDS unchanged"),passed=TRUE),
  file.path(out,"negative_and_invariance_checks.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
cat("PASS: all RR FC merge regression checks\n")
