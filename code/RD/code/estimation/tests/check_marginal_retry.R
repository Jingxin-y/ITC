## Run from the project root. Controlled first-attempt failures test the
## retry/RNG/bootstrap contract, followed by real optim() on the same data.
## Additional rescue is disabled here to isolate the original retry contract.
source('load_anchored.R',encoding='UTF-8'); source_anchored()
check <- function(ok,label) {stopifnot(ok);cat('PASS:',label,'\n')}
dat <- expand.grid(id=1:40,x=0:1,z=0:1)
dat$y <- as.integer(dat$id <= 40*c(.2,.3,.3,.6)[1+dat$x+2*dat$z])
va <- vb <- cbind(1,dat$z); w <- ifelse(dat$z==0,1.4,.6)
pars <- c(log(1.5),log(2/1.5),qlogis(.2)+qlogis(.3),log(6))
args <- list(y=dat$y,x=dat$x,va=va,vb=vb,pars=pars,
 weight=rep(1,nrow(dat)),standardization.weight=w,
 compute.interval=FALSE,compute.p.value=FALSE,psi.grid=.1,
 nsim.refine=4L,seed=801L,time.limit.sec=60,
 optim.multistart=FALSE,optim.fallback=FALSE)
baseline <- do.call(exact_marginal_target,args)
check(baseline$diagnostics$counts$retried.profile.fits==0L,
 'regular fits do not perform additional optimizer calls')
## Bind a wrapper only in a private copy of the solver environment.
solver.env <- new.env(parent=globalenv())
sys.source('estimation/helpers/CI_marginal_target.R',solver.env)
real.optim <- stats::optim; calls <- list(); pending <- NULL
solver.env$optim <- function(par,fn,gr,method,lower,upper,control) {
 val <- fn(par); grad <- gr(par)
 calls[[length(calls)+1L]] <<- list(value=val,gradient=grad,control=control)
 if (is.null(pending)) {
  raw <- real.optim(par,fn,gr,method=method,lower=lower,upper=upper,control=control)
  pending <<- list(par=raw$par,value=fn(raw$par),gradient=gr(raw$par))
  raw$convergence <- 1L
  return(raw)
 }
 ## The next attempt must use the same objective/response/psi at the last
 ## parameter vector, with more iterations and a tighter factr.
 check(identical(par,pending$par) && identical(val,pending$value) &&
  identical(grad,pending$gradient) && control$maxit==1000L &&
  control$factr < calls[[length(calls)-1L]]$control$factr &&
  control$pgtol < calls[[length(calls)-1L]]$control$pgtol,
  'retry preserves data/psi/objective and continues the previous parameters')
 pending <<- NULL
 real.optim(par,fn,gr,method=method,lower=lower,upper=upper,control=control)
}
set.seed(90); rng <- .Random.seed
rescued <- do.call(solver.env$exact_marginal_target,args)
check(identical(rng,.Random.seed),'retry preserves caller RNG state')
counts <- rescued$diagnostics$counts
check(counts$bootstrap.draws==4L && counts$failed.bootstrap.draws==0L &&
 counts$recovered.profile.fits==counts$profile.fits &&
 counts$optimizer.attempts==2L*counts$profile.fits,
 'each failed fit is retried; generated bootstrap sample count stays four')
a <- baseline$diagnostics$evaluations; b <- rescued$diagnostics$evaluations
check(identical(a$acceptability,b$acceptability) &&
 max(abs(a$statistic-b$statistic))<1e-6 &&
 abs(baseline$estimate-rescued$estimate)<1e-12,
 'retry does not change the observed target or the bootstrap sampling distribution')
check(length(rescued$diagnostics$retried.fits)>0L &&
 all(vapply(rescued$diagnostics$retried.fits,function(z)
  length(z$response)==nrow(dat) && length(z$fit$attempts)==2L,logical(1))),
 'retried responses and attempt diagnostics are retained')
## Persistent failure remains a failure; none of its bootstrap draws disappear.
solver.env$optim <- function(...) {raw<-real.optim(...);raw$convergence<-1L;raw}
failed <- do.call(solver.env$exact_marginal_target,args)
check(failed$diagnostics$counts$failed.profile.fits>0L &&
 failed$diagnostics$anchor$optimizer.attempts==3L &&
 !failed$diagnostics$anchor$converged &&
 length(failed$diagnostics$failed.fits[[1]]$response)==nrow(dat),
 'exhausted retries remain failed and retain the response')
## A single unrecovered bootstrap fit invalidates the candidate; all B draws
## must still be counted, rather than redrawn or excluded from the denominator.
ncall <- 0L
solver.env$optim <- function(...) {
 ncall <<- ncall+1L; raw <- real.optim(...)
 if (ncall %in% 3:5) raw$convergence <- 1L
 raw
}
unrecovered <- do.call(solver.env$exact_marginal_target,args)
check(unrecovered$diagnostics$counts$bootstrap.draws==4L &&
 unrecovered$diagnostics$counts$failed.bootstrap.draws==1L &&
 unrecovered$diagnostics$evaluations$status=='bootstrap_fit_failed' &&
 is.na(unrecovered$diagnostics$evaluations$acceptability),
 'an unrecovered bootstrap draw stays in the failure count and invalidates acceptance')
## Ordinary numerical errors also retry from a valid evaluation without RNG use.
ncall <- 0L
solver.env$optim <- function(par,fn,gr,...) {
 ncall <<- ncall+1L
 if (ncall==1L) {fn(par);stop('transient numerical error')}
 real.optim(par,fn,gr,...)
}
errored <- do.call(solver.env$exact_marginal_target,args)
check(errored$diagnostics$anchor$converged &&
 errored$diagnostics$anchor$attempts[[1]]$optimizer.code==99L &&
 errored$diagnostics$counts$failed.bootstrap.draws==0L,
 'a transient optimizer error resumes from the last valid evaluation')
## Timeout during a retry must escape the numerical error handler.
ncall <- 0L
solver.env$optim <- function(...) {
 ncall <<- ncall+1L
 if (ncall==2L) stop(structure(list(message='test timeout',call=NULL),
  class=c('marginal_exact_timeout','error','condition')))
 raw<-real.optim(...);raw$convergence<-1L;raw
}
timeout <- tryCatch(do.call(solver.env$exact_marginal_target,args),error=identity)
check(inherits(timeout,'marginal_exact_timeout'),'timeout during retry propagates')
for (bad in list(-1L,1.5,11L,NA_real_)) {
 args$optim.retries <- bad
 check(inherits(try(do.call(exact_marginal_target,args),silent=TRUE),'try-error'),
  'invalid retry count rejected')
}
validate_marginal_bc_control(list(optim.retries=2L,optim.retry.maxit=1000L))
cat('PASS: main-pipeline adapter accepts retry controls\n')
