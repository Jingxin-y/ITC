## Run from the project root. Verify actual profile/bootstrap likelihoods as
## well as the standalone entry points; this is not a coverage simulation.
source("load_anchored.R",encoding="UTF-8")
source_anchored()
check <- function(ok,msg) {stopifnot(ok);cat("PASS:",msg,"\n")}
near <- function(a,b,tol=1e-8)
  isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE))

dat <- expand.grid(id=1:50,x=0:1,z1=0:1)
dat$y <- as.integer(dat$id <= 50*c(.2,.3,.3,.6)[1+dat$x+2*dat$z1])
ipd <- list(data=dat); ad <- get.aggre(ipd); ad$anchored.mean.x[] <- .4
base <- suppressWarnings(get.compare.anchored(ipd,ad,run_exact=FALSE,
  effect_modifiers="z1",model_covariates="z1",weight_covariates="z1",return_full=TRUE))
ones <- rep(1,nrow(dat)); w <- base$weights$weight
check(base$weights$status=="available" && sd(w)>.1 &&
  all(base$model_weights$weight==1),"fixture has unit fits and nonuniform MAIC target weights")

## Independent probability/likelihood calculations inspect every actual fit,
## including both fits on each bootstrap sample.
f <- base$fits$brm
psi <- marginal_logrr_target(f$point.est,f$va,f$vb,w)
check(abs(psi-marginal_logrr_target(f$point.est,f$va,f$vb,ones))>1e-3,
  "MAIC standardization differs from the unit-weight population target")
real_problem <- .marginal_exact_problem
checked <- list()
.marginal_exact_problem <- function(...) {
  problem <- real_problem(...)
  original_profile <- problem$profile
  problem$profile <- function(psi,y.local=dat$y) {
    fit <- original_profile(psi,y.local)
    check(fit$converged,"observed/bootstrap constrained fit converges")
    probs <- getProbRR(drop(f$va %*% fit$alpha),drop(f$vb %*% fit$beta))
    risk <- probs[cbind(seq_len(nrow(dat)),dat$x+1L)]
    independent <- -sum(dbinom(y.local,1,risk,log=TRUE))
    check(near(fit$value,independent,tol=1e-7),
      "constrained likelihood equals the independent unit-weight Bernoulli likelihood")
    check(abs(marginal_logrr_target(fit$par,f$va,f$vb,w)-psi)<1e-7,
      "each observed/bootstrap fit satisfies the MAIC marginal constraint")
    checked[[length(checked)+1L]] <<- list(value=fit$value,
      weighted.value=-sum(w/mean(w)*dbinom(y.local,1,risk,log=TRUE)))
    fit
  }
  problem
}
ci <- tryCatch(exact_marginal_target(y=dat$y,x=dat$x,va=f$va,vb=f$vb,
  pars=f$point.est,weight=ones,standardization.weight=w,
  compute.interval=FALSE,compute.p.value=FALSE,psi.grid=psi+c(-.2,.2),
  nsim.refine=6L,seed=901L,time.limit.sec=60),
  finally={.marginal_exact_problem <- real_problem})
check(ci$diagnostics$counts$bootstrap.draws==12L &&
  ci$diagnostics$counts$failed.bootstrap.draws==0L && length(checked)==27L,
  "both bootstrap constraints were checked on all twelve generated samples")
check(any(vapply(checked,function(z)abs(z$value-z$weighted.value)>1e-3,logical(1))),
  "the checked likelihoods are distinguishable from MAIC-weighted likelihoods")

## Capture the public solver arguments from a standalone BRM/FC call. The
## caller's MAIC weight is preserved for the target before unit fitting begins.
real_solver <- exact_marginal_target; calls <- list()
exact_marginal_target <- function(...) {
  a <- list(...); calls[[length(calls)+1L]] <<- a
  point <- marginal_logrr_target(a$pars,a$va,a$vb,a$standardization.weight)
  list(estimate=point,low=point-.5,up=point+.5,p=.5,status="available",
    interval.available=TRUE,p.status="available",diagnostics=list())
}
standalone <- tryCatch(suppressWarnings(get.brm.bc(ipd,ad,weight=w,
  gamma=base$weights$ga,brm.est=f,bayes.est=make_na_est(),
  effect_modifiers="z1",model_covariates="z1",run_exact=TRUE)),
  finally={exact_marginal_target <- real_solver})
check(length(calls)==2L && all(vapply(calls,function(a)
  all(a$weight==1) && near(a$standardization.weight,w),logical(1))),
  "standalone BRM and FC intervals default to unit likelihood and MAIC target weights")
check(near(standalone$firth$point.est,base$fits[["brm-FC"]]$point.est),
  "standalone FC point fit uses the same unit-weight model as the main pipeline")

## Exercise the real FC example, whose former CI likelihood used MAIC weights.
source("estimation/example_marginal_exact.R",encoding="UTF-8")
example <- suppressWarnings(example_marginal_exact(method="brm-FC",n=200L,seed=1L,
  ci_control=list(compute.interval=FALSE,compute.p.value=FALSE)))
check(all(example$model_weights$weight==1) && sd(example$weights$weight)>.1 &&
  near(example$ci$estimate,marginal_logrr_target(example$fit$point.est,
    example$fit$va,example$fit$vb,example$weights$weight)),
  "FC example keeps unit fitting weights separate from nonuniform MAIC target weights")
blocked <- try(validate_marginal_bc_control(list(weight=w)),silent=TRUE)
check(inherits(blocked,"try-error"),"example tuning cannot overwrite likelihood or target weights")
cat("Unit likelihood / MAIC target checks passed; confidence-interval coverage was not evaluated.\n")
