## Regression fixtures distilled from failed simulated common/rare records.
## Run from the code root: Rscript estimation/tests/check_brm_fc_solver.R
source("load_anchored.R",encoding="UTF-8"); source_anchored()
fixtures <- list(
  common=list(h=c(-1.31886405305408472,-.06280305014543265,
                   1.19325795276321922,2.44931895567187130),
    n0=c(13,21,9,7),s0=c(7,11,4,1),n1=c(9,22,19,0),s1=c(3,12,10,0)),
  rare=list(h=c(-1.0638689687335476,.2828006119418291,
                 1.6294701926172057,2.9761397732925827),
    n0=c(7,16,18,9),s0=c(0,0,0,0),n1=c(1,16,27,6),s1=c(0,0,2,3)),
  rare_multistart_164=list(h=c(-1.06386005784919,.24954742097697,
    1.56295489980313,2.87636237862929),n0=c(7,14,21,8),s0=c(1,2,0,1),
    n1=c(4,26,14,6),s1=c(1,0,0,0)),
  rare_multistart_155=list(h=c(-1.08667979784149,.0817931030633378,
    1.25026600396817,2.41873890487299),n0=c(12,25,9,4),s0=c(1,0,0,0),
    n1=c(12,19,16,3),s1=c(0,2,0,0)),
  rare_multistart_409=list(h=c(-1.19314767262531,.261910464722629,
    1.71696860207057,3.1720267394185),n0=c(11,13,24,2),s0=c(0,2,0,0),
    n1=c(7,24,17,2),s1=c(1,0,0,0)))
for (name in names(fixtures)) {
  fixture <- fixtures[[name]]
  blocks <- list()
  for (cell in 1:4) for (arm in 0:1) {
    n <- fixture[[paste0("n",arm)]][cell]
    s <- fixture[[paste0("s",arm)]][cell]
    if(n) blocks[[length(blocks)+1L]] <- data.frame(cell=cell,arm=arm,
      h=fixture$h[cell],y=c(rep(1,s),rep(0,n-s)))
  }
  dat <- do.call(rbind,blocks)
  va <- cbind(1,dat$h)
  vb <- model.matrix(~factor(cell,levels=1:4),dat)
  args <- list(param="RR",y=dat$y,x=dat$arm,va=va,vb=vb,weight=rep(1,nrow(dat)),
    alpha.start=c(0,0),beta.start=rep(0,4),max.step=300L,thres=1e-6,pa=2,pb=4)
  fit <- do.call(max.likelihood.firth.rr,args)
  stopifnot(isTRUE(fit$convergence),fit$step<=300L)
  if(grepl("multistart",name)) {
    stopifnot(fit$diagnostics$restarts>1L)
    history <- fit$diagnostics$restart_history
    stopifnot(sum(vapply(history,function(h)h$iterations,integer(1)))==fit$step)
    extra <- Filter(function(h)!h$start%in%c("supplied","pooled_op"),history)
    stopifnot(all(vapply(extra,function(h)h$iterations<=50L,logical(1))))
  }
  terms <- firth_rr_derivatives(fit$par[1:2],fit$par[3:6],dat$arm,va,vb)
  check <- firth_score_adjustment(dat$y,va,vb,args$weight,terms)
  stopifnot(max(abs(check$adjusted.score/fit$diagnostics$score_scale))<=1e-6,
    all(is.finite(check$covariance)),all(diag(check$covariance)>0),
    min(eigen(check$covariance,symmetric=TRUE,only.values=TRUE)$values)>0)
  args$weight <- 19*args$weight
  scaled <- do.call(max.likelihood.firth.rr,args)
  stopifnot(isTRUE(scaled$convergence),max(abs(fit$par-scaled$par))<1e-9)
  args$weight <- args$weight/19
  wrapped <- do.call(MLEst,c(args,list(group="AD",gamma.weight=0,target_mean=0,method="firth")))
  stopifnot(isTRUE(wrapped$converged),all(is.finite(wrapped$se.est)),
            max(abs(wrapped$cov-check$covariance))<1e-9)
  args$max.step <- 1L
  short <- do.call(max.likelihood.firth.rr,args)
  stopifnot(!short$convergence,short$step==1L,short$diagnostics$restarts==0L)
  cat("PASS:",name,"difficult FC root, independent score/covariance, weight scaling and strict budget\n")
}
