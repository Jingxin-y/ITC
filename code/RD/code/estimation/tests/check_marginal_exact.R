## Run from anchored_v1: Rscript estimation/tests/check_marginal_exact.R
## Optional args: project directory, alternate candidate CI_marginal_target.R.
args <- commandArgs(trailingOnly=TRUE)
project <- if(length(args)) args[1] else "."
source(file.path(project,"load_anchored.R"),encoding="UTF-8")
source_anchored(project)
source(if(length(args)>1L) args[2] else
  file.path(project,"estimation/helpers/CI_marginal_target.R"),encoding="UTF-8")
check <- function(ok,message) { stopifnot(ok); cat("PASS:",message,"\n") }
expect_error <- function(expr) inherits(try(force(expr),silent=TRUE),"try-error")

## Stable risk conversion agrees with the existing model in its interior, and
## preserves BOTH defining equations outside the legacy +/-12 approximation.
t <- c(-20,-4,-1,-1e-8,0,1e-8,.3,4,20)
s <- c(-10,2,-2,1,0,-1,.8,-2,-10)
z <- .marginal_exact_risks(t,s)
check(max(abs(z$logp1-z$logp0-t))<1e-12,"risk conversion preserves log RR")
check(max(abs(qlogis(z$p0)+qlogis(z$p1)-s))<1e-5,"risk conversion preserves log OP")
keep <- abs(t)<5
check(max(abs(cbind(z$p0,z$p1)[keep,]-getProbRR(t[keep],s[keep])))<1e-8,
      "stable conversion agrees with getProbRR in its regular range")

dat <- expand.grid(id=1:50,x=0:1,z=0:1)
dat$y <- as.integer(dat$id <= 50*c(.2,.3,.3,.6)[1+dat$x+2*dat$z])
va <- vb <- cbind(1,dat$z)
pars <- c(log(1.5),log(2/1.5),qlogis(.2)+qlogis(.3),log(6))
w <- ifelse(dat$z==0,1.2,.8)
tg <- marginal_logrr_target(pars,va,vb,w,TRUE)
check(abs(tg$value-log(1.75))<1e-12,"marginal target equals the independent risk calculation")
check(abs(tg$value-sum(c(1,.4)*pars[1:2]))>.01,"target is not the average conditional log RR")
numgrad <- sapply(seq_along(pars),function(j) {
  h<-1e-5;p<-m<-pars;p[j]<-p[j]+h;m[j]<-m[j]-h
  (marginal_logrr_target(p,va,vb,w)-marginal_logrr_target(m,va,vb,w))/(2*h)
})
check(max(abs(tg$gradient-numgrad))<1e-8,"full alpha/beta target gradient matches finite differences")
check(abs(marginal_logrr_target(pars,cbind(1,0:1),cbind(1,0:1),c(.6,.4))-tg$value)<1e-12,
      "separate standardization design gives the same population target")

problem <- .marginal_exact_problem(dat$y,dat$x,va,vb,w,pars,va,vb,w,1L,
  8,10,200,1e-10,1e-4,1e-10,200,1e-7,1e-12,function(stage)NULL)
for (psi in c(-.4,tg$value,1.2)) {
  fit <- problem$profile(psi)
  check(fit$converged && abs(marginal_logrr_target(fit$par,va,vb,w)-psi)<1e-7,
        paste("nonlinear constraint after nuisance optimization at",round(psi,4)))
  u <- problem$start+c(.05,-.04,.03)
  a <- problem$evaluate(u,psi,dat$y)
  g <- sapply(seq_along(u),function(j) {
    h<-1e-5;p<-m<-u;p[j]<-p[j]+h;m[j]<-m[j]-h
    (problem$evaluate(p,psi,dat$y)$value-problem$evaluate(m,psi,dat$y)$value)/(2*h)
  })
  check(max(abs(a$gradient-g))<2e-5,"implicit constrained gradient matches numerical likelihood derivative")
}

## Pure marginal limit: constraint must be alpha0=psi, with OP still profiled.
v0 <- matrix(1,nrow(dat),1)
p0 <- c(.2,-1)
pr0 <- .marginal_exact_problem(dat$y,dat$x,v0,v0,w,p0,v0,v0,w,1L,
  8,10,200,1e-10,1e-4,1e-10,200,1e-7,1e-12,function(stage)NULL)
for (psi in c(-.2,.5,1)) {
  f <- pr0$profile(psi)
  independent <- optimize(function(beta) {
    probs<-getProbRR(rep(psi,nrow(dat)),rep(beta,nrow(dat)))
    p<-probs[cbind(seq_len(nrow(dat)),dat$x+1L)]
    -sum(w*(dat$y*log(p)+(1-dat$y)*log1p(-p)))
  },c(-8,8),tol=1e-10)
  check(f$converged && abs(f$alpha-psi)<1e-12 && abs(f$value-independent$objective)<1e-6,
        "no-EM constrained likelihood equals independent coefficient profiling")
}

base <- list(y=dat$y,x=dat$x,va=va,vb=vb,weight=w,pars=pars,psi.se=.2,
  nsim.coarse=30L,nsim.confirm=40L,nsim.refine=40L,nsim.p.value=40L,
  grid.mult=25,thres.dicho=.08,bisection.max.step=10L,seed=902L,
  psi.range=c(-.3,1.3),time.limit.sec=180)
set.seed(888); before <- .Random.seed
ci <- do.call(exact_marginal_target,base)
check(identical(before,.Random.seed),"explicit seed does not alter the caller RNG state")
check(ci$interval.available && ci$low<=ci$estimate && ci$up>=ci$estimate && is.finite(ci$p),
      "complete bootstrap inversion yields a marginal interval and p-value")
check(ci$diagnostics$counts$failed.bootstrap.draws==0L,"no bootstrap fits silently excluded")
check(max(vapply(ci$diagnostics$candidates,function(f)abs(f$constraint.error),numeric(1)))<1e-7,
      "all observed candidates satisfy their marginal constraints")
scaled <- base;scaled$weight<-w*100;scaled$standardization.weight<-w*17
ci2 <- do.call(exact_marginal_target,scaled)
check(max(abs(unlist(ci[c('estimate','low','up','p')])-unlist(ci2[c('estimate','low','up','p')])))<1e-8,
      "reproducible inversion is invariant to scaling either set of weights")

bad <- base;bad$pars<-c(.1,.1,-.3,.1);bad$optim.maxit<-1L
bad$optim.retries<-0L
bad$optim.multistart<-FALSE;bad$optim.fallback<-FALSE # isolate an exhausted original path
bad$compute.p.value<-FALSE
failure <- do.call(exact_marginal_target,bad)
check(!failure$interval.available && is.na(failure$low) && length(failure$diagnostics$failed.fits)>0,
      "nonconverged observed fits return failure plus saved last parameters")
edge <- base;edge$psi.range<-tg$value+c(-1,1)*1e-5;edge$compute.p.value<-FALSE
trunc <- do.call(exact_marginal_target,edge)
check(trunc$interval.available && trunc$status=='boundary_truncated',
      "accepted numerical boundaries are retained and explicitly labelled")
timeout <- base;timeout$time.limit.sec<-1e-12
e <- tryCatch(do.call(exact_marginal_target,timeout),error=function(e)e)
check(inherits(e,'marginal_exact_timeout'),"timeout propagates instead of becoming a numerical fit")
invalid <- base;invalid$va[,1]<-2
check(expect_error(do.call(exact_marginal_target,invalid)),"invalid intercept design is rejected")
invalid<-base;invalid$standardization.weight<-rep(0,nrow(dat))
check(expect_error(do.call(exact_marginal_target,invalid)),"zero standardization mass is rejected")
cat('SMOKE RESULT:',ci$estimate,ci$low,ci$up,ci$p,'\n')
cat('Small bootstrap budgets above verify implementation, NOT 95% coverage.\n')
