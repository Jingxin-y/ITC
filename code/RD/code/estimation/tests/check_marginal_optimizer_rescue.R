## Run from the RR project root; optional argument is the pre-update solver.
## Small fixed problems check numerical behavior, not Monte Carlo coverage.
args <- commandArgs(trailingOnly=TRUE)
source("load_anchored.R",encoding="UTF-8"); source_anchored()
solver_file <- normalizePath("estimation/helpers/CI_marginal_target.R",winslash="/")
check <- function(ok,label) {stopifnot(isTRUE(ok));cat("PASS:",label,"\n")}
private <- function(file=solver_file) {e<-new.env(parent=globalenv());sys.source(file,e);e}
real.optim <- stats::optim; real.nlminb <- stats::nlminb
dat <- expand.grid(id=1:40,x=0:1,z=0:1)
dat$y <- as.integer(dat$id <= 40*c(.2,.3,.3,.6)[1+dat$x+2*dat$z])
y <- dat$y; x <- dat$x; va <- vb <- cbind(1,dat$z)
w <- ifelse(dat$z==0,1.4,.6)
pars <- c(log(1.5),log(2/1.5),qlogis(.2)+qlogis(.3),log(6))
problem <- function(e=.GlobalEnv,...) e$.marginal_exact_problem(
  y,x,va,vb,rep(1,length(y)),pars,va,vb,w,1L,15,15,200L,1e-10,1e-4,
  1e-10,200L,1e-7,1e-12,function(stage)NULL,...)

base <- problem(optim.multistart=FALSE,optim.fallback=FALSE)$profile(.1)
regular <- problem()$profile(.1)
check(base$converged && regular$converged && identical(base$par,regular$par) &&
  identical(base$value,regular$value) && regular$optimizer.attempts==1L &&
  regular$multistart.attempts==0L && regular$fallback.attempts==0L,
  "original successful optimizer path and solution are unchanged")
if(length(args)) {
  old_env <- private(normalizePath(args[1],winslash="/"))
  old <- problem(old_env)$profile(.1)
  check(old$converged && identical(old$par,regular$par) && identical(old$value,regular$value),
    "pre-update solver and new successful profile are exactly identical")
  assignment <- function(f,key) {
    expressions <- as.list(body(f))[-1L]
    selected <- Filter(function(z)is.call(z) && identical(z[[1L]],as.name("<-")) &&
      identical(z[[2L]],as.name(key)),expressions)
    stopifnot(length(selected)==1L);selected[[1L]]
  }
  check(all(vapply(c("target","unpack","evaluate"),function(key)
      identical(assignment(old_env$.marginal_exact_problem,key),
                assignment(.marginal_exact_problem,key)),logical(1))),
    "target, marginal constraint, likelihood and analytic gradient definitions are unchanged")
  check(all(vapply(c("evaluate.acceptability","endpoint"),function(key)
      identical(assignment(old_env$exact_marginal_target,key),
                assignment(exact_marginal_target,key)),logical(1))),
    "bootstrap sampling/statistic and interval endpoint search definitions are unchanged")
}

## Only the original attempt and its two resumes report failure.
e <- private(); ncall<-0L; starts<-list(); p<-NULL
e$optim <- function(par,fn,gr,...) {
  ncall<<-ncall+1L; starts[[ncall]]<<-par
  exact <- p$evaluate(par,.1,y)
  stopifnot(abs(fn(par)-exact$value/length(y))<1e-12,
    max(abs(gr(par)-exact$gradient/length(y)))<1e-12)
  raw <- real.optim(par,fn,gr,...)
  if(ncall<=3L) raw$convergence<-1L
  raw
}
p <- problem(e); set.seed(611L); rng <- .Random.seed
z <- p$profile(.1)
check(z$converged && z$retry.attempts==3L && z$multistart.attempts>0L &&
  z$fallback.attempts==0L && any(vapply(starts[-seq_len(3L)],
    function(v)max(abs(v-starts[[1L]]))>1e-6,logical(1))),
  "distinct deterministic starts run only after all three original attempts fail")
check(identical(rng,.Random.seed) && abs(z$value-base$value)<1e-6 &&
  abs(z$constraint.error)<=1e-7 && z$max.projected.score<=1e-4 &&
  all(vapply(starts,function(v)all(abs(v)<=15),logical(1))),
  "multistart preserves RNG, response, likelihood, bounds and fit tolerances")

## Both kinds of neighboring successful fits are available for failed profiles.
e <- private(); p <- problem(e); invisible(p$profile(.15))
e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
z <- p$profile(.1)
check(z$converged && z$fallback.attempts>0L && z$optimizer.method=="nlminb" &&
  any(vapply(z$attempts,function(a)a$start.label=="same_response_neighbor",logical(1))),
  "same-response neighboring solution supplies a start and nlminb recovers failure")

## A zero status code and low likelihood alone must not pass the gradient test.
e <- private()
e$optim <- function(par,fn,gr,...) list(par=rep(0,length(par)),convergence=0L,counts=c(0L,0L))
e$nlminb <- function(start,objective,gradient,...) list(par=rep(0,length(start)),
  convergence=0L,evaluations=c(0L,0L))
z <- problem(e)$profile(.1)
check(!z$converged && z$max.projected.score>1e-4,
  "all optimizers remain subject to the original independent acceptance checks")

arguments <- list(y=y,x=x,va=va,vb=vb,pars=pars,weight=rep(1,length(y)),
  standardization.weight=w,compute.interval=FALSE,compute.p.value=FALSE,
  psi.grid=.1,nsim.refine=4L,seed=731L,time.limit.sec=60)
e <- private(); failed.response <- NULL
is.failed <- function(fn) {
  yy <- base::get("y.local",envir=environment(fn),inherits=TRUE)
  if(!identical(yy,y) && is.null(failed.response)) failed.response<<-yy
  !is.null(failed.response) && identical(yy,failed.response)
}
e$optim <- function(par,fn,gr,...) {
  bad<-is.failed(fn);raw<-real.optim(par,fn,gr,...)
  if(bad)raw$convergence<-1L;raw
}
e$nlminb <- function(start,objective,gradient,...) {
  bad<-is.failed(objective);raw<-real.nlminb(start,objective,gradient=gradient,...)
  if(bad)raw$convergence<-1L;raw
}
z <- do.call(e$exact_marginal_target,arguments)
check(z$diagnostics$counts$bootstrap.draws==4L &&
  z$diagnostics$counts$failed.bootstrap.draws==1L &&
  z$diagnostics$evaluations$status=="bootstrap_fit_failed" &&
  is.na(z$diagnostics$evaluations$acceptability),
  "an exhausted bootstrap response remains failed and is never dropped or redrawn")
check(length(z$diagnostics$failed.fits)>0L &&
  all(vapply(z$diagnostics$failed.fits,function(a)length(a$response)==length(y),logical(1))) &&
  z$diagnostics$counts$multistart.attempts>0L && z$diagnostics$counts$fallback.attempts>0L,
  "failed responses, stage, attempts and recovery counters are retained")

## Exercise the complete candidate/anchor/bootstrap adapter through nlminb.
e <- private();e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
set.seed(912L);rng<-.Random.seed
a <- do.call(exact_marginal_target,arguments)
b <- do.call(e$exact_marginal_target,arguments)
check(identical(rng,.Random.seed) && b$diagnostics$counts$bootstrap.draws==4L &&
  b$diagnostics$counts$failed.bootstrap.draws==0L &&
  b$diagnostics$counts$recovered.by.fallback>0L &&
  abs(a$estimate-b$estimate)<1e-12 &&
  max(abs(a$diagnostics$evaluations$statistic-b$diagnostics$evaluations$statistic))<1e-6 &&
  identical(a$diagnostics$evaluations$acceptability,b$diagnostics$evaluations$acceptability),
  "fallback retains the original draws, bootstrap statistic, acceptability and input estimate")

## A success on the first resume also bypasses all new rescue calls.
e <- private(); calls<-0L
e$optim <- function(...) {calls<<-calls+1L;raw<-real.optim(...);if(calls==1L)raw$convergence<-1L;raw}
e$nlminb <- function(...) stop("fallback must not run after an original resume succeeds")
z <- problem(e)$profile(.1)
check(z$converged && z$retry.attempts==2L && z$multistart.attempts==0L &&
  z$fallback.attempts==0L,"success on an original resume bypasses additional rescue")

e <- private();e$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
e$nlminb <- function(...) stop(structure(list(message="test timeout",call=NULL),
  class=c("marginal_exact_timeout","error","condition")))
z <- tryCatch(problem(e)$profile(.1),error=identity)
check(inherits(z,"marginal_exact_timeout"),"fallback timeout propagates")

for(ctrl in list(list(optim.multistart=NA),list(optim.multistart=1),list(optim.fallback="yes"),
  list(optim.multistart.max=0L),list(optim.multistart.max=9L),list(optim.multistart.max=1.5))) {
  z<-tryCatch(do.call(exact_marginal_target,c(arguments,ctrl)),error=identity)
  check(inherits(z,"error"),paste("invalid rescue control rejected:",names(ctrl)))
}
validate_marginal_bc_control(list(optim.multistart=TRUE,optim.multistart.max=3L,optim.fallback=TRUE))
check(TRUE,"main pipeline accepts the new optional rescue controls")

## Whole-interval regression: use identical generated responses and unchanged searches.
if(length(args)) {
  old <- private(normalizePath(args[1],winslash="/"))
  ci_args <- list(y=y,x=x,va=va,vb=vb,pars=pars,weight=w,
    psi.se=.2,psi.range=c(-.3,1.3),nsim.coarse=12L,nsim.confirm=12L,
    nsim.refine=12L,nsim.p.value=12L,grid.mult=25,thres.dicho=.08,
    bisection.max.step=10L,seed=902L,time.limit.sec=180)
  before <- do.call(old$exact_marginal_target,ci_args)
  after <- do.call(exact_marginal_target,ci_args)
  check(before$interval.available && after$interval.available &&
    identical(before[c("estimate","low","up","p","status","p.status")],
              after[c("estimate","low","up","p","status","p.status")]) &&
    identical(before$diagnostics$evaluations,after$diagnostics$evaluations) &&
    after$diagnostics$counts$multistart.attempts==0L && after$diagnostics$counts$fallback.attempts==0L,
    "pre-update successful complete interval, P value and bootstrap search are exactly unchanged")
}
cat("PASS: RR OPTIMIZER RESCUE CHECKS COMPLETE\n")
