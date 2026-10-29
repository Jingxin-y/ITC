## RR and standalone RD must share the FC equation, solver and controls.
## Their effect links and probability derivatives intentionally differ.
args <- commandArgs(trailingOnly=TRUE)
rd_dir <- if(length(args)) args[1L] else "RD"
source(file.path(rd_dir,"load_anchored_rd.R"),encoding="UTF-8")
source_anchored_rd(".",rd_dir)
check <- function(ok,label) {stopifnot(isTRUE(ok));cat("PASS:",label,"\n")}
private <- .anchored_rd_solver_env
for(key in c("get.estimate","MLEst","max.likelihood.firth","max.likelihood.firth.rd")) {
  global <- get(key,envir=.GlobalEnv)
  local <- get(key,envir=private,inherits=FALSE)
  check(identical(body(global),body(local)) && identical(formals(global),formals(local)),
    paste("shared function body and arguments:",key))
}
check(identical(get("firth_score_adjustment",envir=private,inherits=TRUE),firth_score_adjustment),
  "RR and RD resolve the same fixed-analytic-weight FC score adjustment")
controls <- anchored_fit_control()
check(controls$fc_max_step==300L && controls$fc_thres==1e-6,
  "common default FC total budget is 300 and scaled score tolerance is 1e-6")
dat <- expand.grid(id=1:40,x=0:1,z1=0:1)
dat$y <- as.integer(dat$id <= 40*c(.2,.3,.3,.6)[1+dat$x+2*dat$z1])
ipd <- list(data=dat)
for(mode in c("unit","analytic")) {
  w <- if(mode=="unit") rep(1,nrow(dat)) else exp(.4*dat$z1)
  for(model in c("RR","RD")) {
    fit <- suppressWarnings(if(model=="RR") get.estimate("AD",ipd,w,
        gamma.ipd=numeric(0),target_mean=numeric(0),method="firth",param="RR",
        effect_modifiers="z1",model_covariates="z1",non_effect_modifiers=character(0),
        weight_covariates=character(0),op_model="marginal") else
      rd_get_estimate("AD",ipd,w,gamma.ipd=numeric(0),target_mean=numeric(0),
        method="firth",param="RD",effect_modifiers="z1",model_covariates="z1",
        non_effect_modifiers=character(0),weight_covariates=character(0),op_model="marginal"))
    check(fit$converged && fit$fc.diagnostics$maximum_iterations==300L &&
      fit$fc.diagnostics$score_tolerance==1e-6 && fit$fc.diagnostics$weights_normalized_to_mean_one,
      paste(model,mode,"pipeline uses the shared FC budget, tolerance and weight normalization"))
    derivative <- if(model=="RR") firth_rr_derivatives else private$firth_rd_derivatives
    nalpha <- ncol(fit$va); nbeta <- ncol(fit$vb)
    call <- list(y=dat$y,x=dat$x,va=fit$va,vb=fit$vb,alpha.start=rep(0,nalpha),
      beta.start=rep(0,nbeta),weight=w,max.step=300L,thres=1e-6,
      pa=nalpha,pb=nbeta,derivatives=derivative)
    direct <- do.call(max.likelihood.firth,call)
    check(direct$convergence && identical(unname(direct$par),unname(fit$point.est)) &&
      identical(direct$diagnostics,fit$fc.diagnostics),
      paste(model,mode,"pipeline exactly reproduces the same core solver with its model derivative"))
    normalized <- w/max(w); normalized <- normalized/mean(normalized)
    terms <- derivative(fit$point.est[seq_len(nalpha)],
      fit$point.est[nalpha+seq_len(nbeta)],dat$x,fit$va,fit$vb)
    equation <- firth_score_adjustment(dat$y,fit$va,fit$vb,normalized,terms)
    check(max(abs(equation$adjusted.score/fit$fc.diagnostics$score_scale))<=1e-6,
      paste(model,mode,"accepted solution independently satisfies the same adjusted-score test"))
    call$weight <- 19*w
    scaled <- do.call(max.likelihood.firth,call)
    check(scaled$convergence && max(abs(scaled$par-direct$par))<1e-8,
      paste(model,mode,"FC root is invariant to global weight rescaling"))
    call$weight <- w;call$max.step<-1L
    short <- do.call(max.likelihood.firth,call)
    check(!short$convergence && short$step==1L,
      paste(model,mode,"unsolved adjusted score is rejected at a strict one-step budget"))
  }
}
check(identical(private$probability_mapping,"stable_rd_op_v1"),
  "standalone RD deliberately rebinds only its stable RD/OP probability mapping")
cat("PASS: RR/RD FIRTH SOLVER CONSISTENCY CHECKS COMPLETE\n")
