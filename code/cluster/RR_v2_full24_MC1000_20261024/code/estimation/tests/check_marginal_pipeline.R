## Run from anchored_v1. Small B verifies integration, not coverage accuracy.
source("load_anchored.R",encoding="UTF-8")
source_anchored()
check <- function(ok,msg) {stopifnot(ok);cat("PASS:",msg,"\n")}
equal <- function(a,b) isTRUE(all.equal(a,b,tolerance=1e-10))
dat <- expand.grid(id=1:50,x=0:1,z1=0:1)
dat$y <- as.integer(dat$id <= 50*c(.2,.3,.3,.6)[1+dat$x+2*dat$z1])
ipd <- list(data=dat); ad <- get.aggre(ipd); ad$anchored.mean.x[] <- .4
args <- list(data.ipd=ipd,data.ad=ad,effect_modifiers="z1",model_covariates="z1",
  weight_covariates="z1",return_full=TRUE,seed=902L)
control <- list(nsim.coarse=6L,nsim.confirm=6L,nsim.refine=6L,nsim.p.value=6L,
  thres.dicho=.25,grid.mult=25,bisection.max.step=10L,psi.range=c(-.6,1.6),time.limit.sec=60)
old_exact <- exact
exact <- function(...) stop("The coefficient exact() must not run in marginal mode.")
set.seed(77); before <- .Random.seed
result <- suppressWarnings(do.call(get.compare.anchored,c(args,list(run_exact=TRUE,bc_control=control))))
check(identical(before,.Random.seed),"seeded marginal BC preserves caller RNG")
rows <- result$marginal_estimates
check(nrow(rows)==3L && all(rows$component=="marginal") && all(rows$parameter=="log_RR") &&
  all(is.na(rows$se)),"BC rows contain marginal targets and do not recycle coefficient SEs")
check(all(subset(rows,method=="brm_b-BC")$status=="not_applicable_with_effect_modifiers"),
  "scalar Bayes limitation remains explicit")
for (j in 1:2) {
  name <- c("brm-BC","brm-FC-BC")[j]; f <- result$fits[[name]]
  direct <- do.call(exact_marginal_target,c(list(y=dat$y,x=dat$x,va=f$va,vb=f$vb,
    pars=f$point.est,weight=result$model_weights$weight,standardization.weight=result$weights$weight,
    seed=902L+1000000L),control))
  actual <- result$bc[[j]]
  check(equal(actual[c("estimate","low","up","p","status")],
             direct[c("estimate","low","up","p","status")]),paste(name,"equals direct exact_marginal_target"))
  row <- subset(rows,method==name)
  check(equal(row$estimate,actual$estimate) && actual$interval.available &&
    equal(c(row$lower,row$upper,row$p.value),c(actual$low,actual$up,actual$p)),
    paste(name,"assembles genuine bootstrap endpoints and p-value"))
}
off <- suppressWarnings(do.call(get.compare.anchored,c(args,list(run_exact=FALSE))))
check(equal(result$fits$brm,off$fits$brm) &&
  equal(result$fits[["brm-FC"]],off$fits[["brm-FC"]]) &&
  all(off$marginal_estimates$interval=="BC_not_run") && all(is.na(off$marginal_estimates$lower)),
  "disabling BC retains fitted models and marginal points without intervals")

## Failure and boundary statuses must propagate independently of the point estimate.
f <- result$fits$brm; w <- result$weights$weight
boundary <- result$bc[[1]]; boundary$status <- "boundary_truncated"
b <- comparison_marginal_rows("brm-BC",f,boundary,w,TRUE)
check(b$interval_status=="BC_boundary_truncated" && is.finite(b$lower),"boundary truncation is preserved")
failed <- result$bc[[1]]; failed$status <- "interval_failed"; failed$interval.available <- FALSE
failed$low <- NA_real_
b <- comparison_marginal_rows("brm-BC",f,failed,w,TRUE)
check(b$status=="available" && is.finite(b$estimate) && is.na(b$lower) && is.na(b$upper),
  "interval failure does not erase the marginal point or keep a half interval")
timeout <- suppressWarnings(marginal_exact_from_fit(f,dat$y,dat$x,w,w,TRUE,902L,
  list(time.limit.sec=1e-12,psi.range=c(-.6,1.6))))
check(timeout$status=="timeout" && is.finite(timeout$estimate) && is.na(timeout$low),
  "timeout retains point and error diagnostics")
invalid <- try(do.call(get.compare.anchored,c(args,list(bc_control=list(weight=rep(1,nrow(dat)))))),silent=TRUE)
check(inherits(invalid,"try-error"),"CI tuning cannot silently replace the fitted data or target weights")

## No EM: scalar matrix output and Bayes covariance with unavailable cross terms.
scalar <- suppressWarnings(get.compare.anchored(ipd,ad,run_exact=FALSE,weight_covariates="z1",
  model_covariates="z1",effect_modifiers=character(0),return_full=FALSE))
check(identical(dim(scalar),c(5L,length(anchored_method_names()))) && nrow(attr(scalar,"marginal_estimates"))==3L &&
  all(is.finite(scalar["point.est",c("brm-BC","brm-FC-BC","brm_b-BC")])),
  "scalar matrix includes marginal BC points and preserves its public shape")
bad_ad <- ad; bad_ad$anchored.mean.x[] <- 2
bad <- suppressWarnings(do.call(get.compare.anchored,modifyList(args,list(data.ad=bad_ad,run_exact=TRUE))))
check(all(subset(bad$marginal_estimates,method!="brm_b-BC")$status=="weighting_failed"),
  "weighting failure preserves marginal schema")

## Exercise run/simulate/summary forwarding without a costly second bootstrap run.
real_marginal <- exact_marginal_target; calls <- 0L
exact_marginal_target <- function(...) {
  a <- list(...); calls <<- calls+1L
  p <- marginal_logrr_target(a$pars,a$va,a$vb,a$standardization.weight)
  list(estimate=p,low=p-.5,up=p+.5,p=.5,status="available",interval.available=TRUE,
       p.status="available",diagnostics=list(marker="marginal forwarding"))
}
batch <- suppressWarnings(simulate_anchored(joint_rr_gop_design(event="common",ess_ratio=.5),
  n_AC=200,n_BC=200,seeds=1L,score_spec=make_brm_score_spec(),
  weight_covariates=paste0("z",1:3),run_exact=TRUE))
check(calls==2L && batch$bc_target=="marginal" &&
  nrow(subset(batch$estimates,component=="marginal"))==3L,
  "run.anchored and simulate_anchored default to marginal exact")
s <- summarize_anchored_simulation(batch)
check(nrow(subset(s$metrics,component=="marginal"))==3L &&
  !any(batch$reference_truth$method %in% c("brm-BC","brm-FC-BC")),
  "summary includes marginal rows without assigning conditional coefficient truths")
truth <- data.frame(method="brm-BC",component="marginal",term="marginal log RR",truth=.2)
s <- summarize_anchored_simulation(batch,truth)
check(subset(s$metrics,method=="brm-BC")$n_coverage==1L,"explicit marginal truth is usable in summaries")
exact_marginal_target <- real_marginal; exact <- old_exact
cat("Small bootstrap budgets validate wiring, not confidence-interval coverage.\n")
