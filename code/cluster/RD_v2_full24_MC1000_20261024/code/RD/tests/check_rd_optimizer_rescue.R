## RD BC numerical rescue checks; no full simulation or coverage claim.
args <- commandArgs(trailingOnly=TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
update_dir <- dirname(dirname(normalizePath(own, winslash="/")))
code_dir <- if (length(args)) args[1] else dirname(update_dir)
rd_dir <- file.path(code_dir,"RD")
solver_file <- if (length(args)>1L) args[2] else file.path(update_dir,"estimation/CI_marginal_rd_target.R")
if (.Platform$OS.type=="windows") Sys.setlocale("LC_CTYPE","English_United States.utf8")
source(file.path(rd_dir,"load_anchored_rd.R"),encoding="UTF-8"); source_anchored_rd(code_dir,rd_dir)
source(solver_file,encoding="UTF-8")
check <- function(ok,msg) {stopifnot(isTRUE(ok));cat("PASS:",msg,"\n")}
private <- function() {e<-new.env(parent=globalenv());sys.source(solver_file,e);e}
real.optim <- stats::optim; real.nlminb <- stats::nlminb
cell <- rep(rep(0:3,each=10L),2L); x <- rep(0:1,each=40L)
va <- vb <- cbind(1,as.numeric(cell==1),as.numeric(cell==2),as.numeric(cell==3))
y <- as.numeric(seq_along(x)%%10L < ifelse(x==0,3L,5L))
w <- exp(.3*cell); pars <- c(atanh(.2),0,0,0,qlogis(.3)+qlogis(.5),0,0,0)
problem <- function(e=.GlobalEnv,...) e$.rd_exact_problem(y,x,va,vb,rep(1,length(y)),pars,
  va,vb,w,1L,15,15,200L,1e-10,1e-4,1e-10,200L,1e-7,1e-12,function(stage)NULL,...)
base <- problem(optim.multistart=FALSE,optim.fallback=FALSE)$profile(.1)
regular <- problem()$profile(.1)
check(base$converged && regular$converged && identical(base$par,regular$par) &&
  identical(base$value,regular$value) && regular$optimizer.attempts==1L,
  "regular fits retain the existing optimizer path and solution")

## Force only the original and resumed starts to report failure.
e <- private(); ncall<-0L; starts<-list(); p<-NULL
e$optim <- function(par,fn,gr,...) {
  ncall<<-ncall+1L; starts[[ncall]]<<-par
  exact <- p$evaluate(par,.1,y)
  stopifnot(abs(fn(par)-exact$value/length(y))<1e-12,
    max(abs(gr(par)-exact$gradient/length(y)))<1e-12)
  raw <- real.optim(par,fn,gr,...)
  if(ncall<=3L)raw$convergence<-1L
  raw
}
p <- problem(e); set.seed(611); rng <- .Random.seed
z <- p$profile(.1)
check(z$converged && z$multistart.attempts>0L && z$fallback.attempts==0L &&
  any(vapply(starts[-seq_len(3L)],function(v)max(abs(v-starts[[1L]]))>1e-6,logical(1))),
  "different deterministic starts recover an exhausted initial retry sequence")
check(identical(rng,.Random.seed) && abs(z$value-base$value)<1e-6 &&
  abs(z$constraint.error)<=1e-7 && z$max.projected.score<=1e-4 &&
  all(vapply(starts,function(v)all(abs(v)<=15),logical(1))),
  "rescue preserves RNG, response, objective, bounds and acceptance tolerances")

## Build a same-response warm start without consuming bootstrap RNG.
e <- private(); p <- problem(e); invisible(p$profile(.15))
e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
z <- p$profile(.1)
check(z$fallback.attempts>0L && z$converged && z$optimizer.method=="nlminb" &&
  any(vapply(z$attempts,function(a)a$start.label=="same_response_neighbor",logical(1))),
  "neighboring constrained fits supply a start and nlminb rescues L-BFGS-B failure")

## Neither a zero optimizer code nor a low objective alone is enough.
e <- private()
e$optim <- function(par,fn,gr,...) list(par=rep(0,length(par)),convergence=0L,counts=c(0L,0L))
e$nlminb <- function(start,objective,gradient,...) list(par=rep(0,length(start)),
  convergence=0L,evaluations=c(0L,0L))
z <- problem(e)$profile(.1)
check(!z$converged && z$max.projected.score>1e-4,
  "all optimizers are rejected when the independent projected score fails")

## A bootstrap response with a permanent failure is never redrawn or dropped.
e <- private(); failed.response <- NULL
is.failed <- function(fn) {
  yy <- base::get("y.local",envir=environment(fn),inherits=TRUE)
  if(!identical(yy,y) && is.null(failed.response)) failed.response<<-yy
  !is.null(failed.response) && identical(yy,failed.response)
}
e$optim <- function(par,fn,gr,...) {
  bad<-is.failed(fn);raw<-real.optim(par,fn,gr,...);if(bad)raw$convergence<-1L;raw
}
e$nlminb <- function(start,objective,gradient,...) {
  bad<-is.failed(objective);raw<-real.nlminb(start,objective,gradient=gradient,...)
  if(bad)raw$convergence<-1L;raw
}
arguments <- list(y=y,x=x,va=va,vb=vb,pars=pars,weight=rep(1,length(y)),
  standardization.weight=w,compute.interval=FALSE,compute.p.value=FALSE,
  psi.grid=.1,nsim.refine=4L,seed=731L,time.limit.sec=60)
z <- do.call(e$exact_marginal_rd_target,arguments)
check(z$diagnostics$counts$bootstrap.draws==4L &&
  z$diagnostics$counts$failed.bootstrap.draws==1L &&
  z$diagnostics$evaluations$status=="bootstrap_fit_failed" &&
  is.na(z$diagnostics$evaluations$acceptability),
  "one exhausted bootstrap draw remains counted and invalidates that candidate")
check(length(z$diagnostics$failed.fits)>0L &&
  all(vapply(z$diagnostics$failed.fits,function(a)length(a$response)==length(y),logical(1))) &&
  z$diagnostics$counts$multistart.attempts>0L && z$diagnostics$counts$fallback.attempts>0L,
  "failed responses, optimizer attempts and rescue counts are retained")

## Test real fallback within the full candidate/anchor/bootstrap adapter.
e <- private();e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
set.seed(912);rng<-.Random.seed
a <- do.call(exact_marginal_rd_target,arguments)
b <- do.call(e$exact_marginal_rd_target,arguments)
check(identical(rng,.Random.seed) && b$diagnostics$counts$bootstrap.draws==4L &&
  b$diagnostics$counts$failed.bootstrap.draws==0L && b$diagnostics$counts$recovered.by.fallback>0L &&
  abs(a$estimate-b$estimate)<1e-12 &&
  max(abs(a$diagnostics$evaluations$statistic-b$diagnostics$evaluations$statistic))<1e-6 &&
  identical(a$diagnostics$evaluations$acceptability,b$diagnostics$evaluations$acceptability),
  "fallback uses the original four draws and preserves the candidate statistic and acceptability")

e <- private();e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
e$nlminb <- function(...) stop(structure(list(message="test timeout",call=NULL),
  class=c("marginal_exact_timeout","error","condition")))
z <- tryCatch(problem(e)$profile(.1),error=identity)
check(inherits(z,"marginal_exact_timeout"),"fallback timeout propagates without being swallowed")
for(ctrl in list(list(optim.multistart=NA),list(optim.multistart=1),list(optim.fallback="yes"),
  list(optim.multistart.max=0L),list(optim.multistart.max=9L),list(optim.multistart.max=1.5))) {
  z<-tryCatch(do.call(exact_marginal_rd_target,c(arguments,ctrl)),error=identity)
  check(inherits(z,"error"),paste("invalid rescue control rejected:",names(ctrl)))
}
validate_rd_bc_control(list(optim.multistart=TRUE,optim.multistart.max=3L,optim.fallback=TRUE))
check(TRUE,"RD comparison adapter accepts all new controls")
cat("PASS: RD OPTIMIZER RESCUE CHECKS COMPLETE\n")
