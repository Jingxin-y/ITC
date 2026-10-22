## Paired numerical profile benchmark, NOT a coverage/BC-return simulation.
args <- commandArgs(trailingOnly=TRUE)
own <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
update_dir <- dirname(dirname(normalizePath(own,winslash="/")))
code_dir <- if(length(args))args[1] else dirname(update_dir)
out <- if(length(args)>1L)args[2] else file.path(update_dir,"validation/paired_profiles")
draws <- if(length(args)>2L)as.integer(args[3]) else 20L
stopifnot(length(draws)==1L,!is.na(draws),draws>=1L)
if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","English_United States.utf8")
rd_dir <- file.path(code_dir,"RD")
source(file.path(rd_dir,"load_anchored_rd.R"),encoding="UTF-8");source_anchored_rd(code_dir,rd_dir)
source(file.path(update_dir,"estimation/CI_marginal_rd_target.R"),encoding="UTF-8")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
new.env <- new.env(parent=globalenv())
sys.source(file.path(update_dir,"estimation/CI_marginal_rd_target.R"),new.env)
old.env <- new.env(parent=globalenv())
old.file <- file.path(update_dir,"validation/CI_marginal_rd_target_before.R")
if(file.exists(old.file)) {
  sys.source(old.file,old.env);old.control<-list()
} else {
  sys.source(file.path(update_dir,"estimation/CI_marginal_rd_target.R"),old.env)
  old.control<-list(optim.multistart=FALSE,optim.fallback=FALSE)
}
cases <- data.frame(event=c("rare","rare","rare","rare","common"),
  strength=c("high","high","low","high","high"),
  weighting=c("unit","unit","unit","maic","unit"),seed=c(19L,29L,41L,19L,19L))
rows<-list();observed<-list();skipped<-list()
field <- function(x,k,default=NA)if(is.null(x[[k]]))default else x[[k]]
make.problem <- function(e,dat,f,wfit,wtarget,control=list())do.call(e$.rd_exact_problem,c(list(
  y=dat$y,x=dat$x,va=f$va,vb=f$vb,weight=wfit,pars=f$point.est,
  va.target=f$va,vb.target=f$vb,standardization.weight=wtarget,
  intercept=1L,alpha.bound=15,beta.bound=15,optim.maxit=200L,
  optim.reltol=1e-10,score.tol=1e-4,root.tol=1e-10,root.maxit=200L,
  constraint.tol=1e-7,eps=1e-12,check.deadline=function(stage)NULL),control))
for(i in seq_len(nrow(cases))) {
  cfg<-cases[i,];design<-webster_rd_scenario_design(cfg$event,.15,"diffuse",cfg$strength)
  fit<-run.anchored.rd(design,100L,seed=cfg$seed,model_weighting=cfg$weighting,run_exact=FALSE)
  dat<-fit$data$data.IPD$data;wt<-fit$weights$weight
  wf<-if(cfg$weighting=="unit")rep(1,nrow(dat)) else wt
  for(method in c("brm","brm-FC")) {
    f<-fit$fits[[method]]
    if(is.null(f)||!isTRUE(f$converged)) {
      skipped[[length(skipped)+1L]]<-data.frame(case=i,method=method,reason="observed point fit unavailable")
      next
    }
    old<-make.problem(old.env,dat,f,wf,wt,old.control);new<-make.problem(new.env,dat,f,wf,wt)
    anchor<-marginal_rd_target(f$point.est,f$va,f$vb,wt)
    targets<-unique(c(anchor,0,.12))
    for(j in seq_along(targets)) {
      psi<-targets[j];a<-old$profile(psi);b<-new$profile(psi)
      observed[[length(observed)+1L]]<-data.frame(case=i,method=method,target=psi,
        old_ok=a$converged,new_ok=b$converged,old_reason=a$reason,new_reason=b$reason)
      ## Exactly the same candidate-specific generating probabilities and yy
      ## are used by BOTH algorithms. Do not discard failed responses.
      generator<-if(isTRUE(a$converged))a else if(isTRUE(b$converged))b else NULL
      if(is.null(generator))next
      pp<-rd_probabilities(drop(f$va%*%generator$alpha),drop(f$vb%*%generator$beta))
      prob<-pmin(pmax(ifelse(dat$x==1,pp$p1,pp$p0),1e-12),1-1e-12)
      set.seed(800000L+i*1000L+j*100L+if(method=="brm")0L else 50L)
      yy<-replicate(draws,rbinom(nrow(dat),1,prob))
      for(k in seq_len(draws))for(target in unique(c(anchor,psi))) {
        t0<-proc.time()[["elapsed"]];a<-old$profile(target,yy[,k]);ta<-proc.time()[["elapsed"]]-t0
        t0<-proc.time()[["elapsed"]];b<-new$profile(target,yy[,k]);tb<-proc.time()[["elapsed"]]-t0
        rows[[length(rows)+1L]]<-data.frame(case=i,event=cfg$event,strength=cfg$strength,
          weighting=cfg$weighting,seed=cfg$seed,method=method,null=psi,target=target,
          draw=k,events=sum(yy[,k]),old_ok=a$converged,new_ok=b$converged,
          old_value=a$value,new_value=b$value,old_seconds=ta,new_seconds=tb,
          old_score=field(a,"max.projected.score"),new_score=field(b,"max.projected.score"),
          new_constraint_error=field(b,"constraint.error"),
          multistart_attempts=field(b,"multistart.attempts",0L),
          fallback_attempts=field(b,"fallback.attempts",0L),
          new_optimizer=field(b,"optimizer.method","L-BFGS-B"),old_reason=a$reason,new_reason=b$reason)
      }
    }
    if(length(rows))write.csv(do.call(rbind,rows),file.path(out,"paired_profiles.csv"),row.names=FALSE)
    cat("COMPLETE case",i,cfg$event,cfg$strength,cfg$weighting,"seed",cfg$seed,method,"\n");flush.console()
  }
}
result<-do.call(rbind,rows)
write.csv(do.call(rbind,observed),file.path(out,"observed_profiles.csv"),row.names=FALSE)
if(length(skipped))write.csv(do.call(rbind,skipped),file.path(out,"skipped_point_fits.csv"),row.names=FALSE)
if(is.null(result)||!nrow(result))stop("No identifiable fitted cases in numerical benchmark.")
result$recovered<-!result$old_ok & result$new_ok;result$lost<-result$old_ok & !result$new_ok
summary<-aggregate(cbind(old_ok,new_ok,recovered,lost,old_seconds,new_seconds)~event+weighting+method,result,mean)
write.csv(summary,file.path(out,"paired_profile_summary.csv"),row.names=FALSE)
print(summary,digits=5)
cat("Profile fits:",nrow(result),"old successes:",sum(result$old_ok),
  "new successes:",sum(result$new_ok),"recovered:",sum(result$recovered),"lost:",sum(result$lost),"\n")
cat("Old/new total solver seconds:",sum(result$old_seconds),sum(result$new_seconds),"\n")
both<-result$old_ok & result$new_ok
cat("Largest NLL improvement on paired successes:",max(result$old_value[both]-result$new_value[both]),"\n")
stopifnot(!any(result$lost),all(abs(result$new_constraint_error[result$new_ok])<=1e-7),
  all(result$new_score[result$new_ok]<=1e-4))
cat("PASS: paired profile checks; no coverage/whole-interval return claim\n")
