## Coefficient-BC compatibility; marginal integration is checked in check_marginal_pipeline.R.
## Run from anchored_v1: Rscript estimation/tests/check_joint_pipeline.R
source("load_anchored.R", encoding="UTF-8")
source_anchored()
close_to <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(a,b,tolerance=tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
design <- joint_rr_gop_design(.4,.2,.3,.5,-3,.2,.4,ess_ratio=.5)
scores <- make_brm_score_spec(c("z1","z4"),c("z2","z3","z5","z6"))
args <- list(bc_target="coefficient",score_standardization="fixed",design=design,n_AC=600,n_BC=800,seed=12,score_spec=scores,
             weight_covariates="h_E",run_exact=FALSE,return_data=TRUE)
result <- suppressWarnings(do.call(run.anchored,args))
set.seed(12)
sim <- generate_anchored_data(600,800,design)
direct <- suppressWarnings(get.compare.anchored(bc_target="coefficient",score_standardization="fixed",sim$data.IPD,sim$data.AD,
  score_spec=scores,weight_covariates="h_E",run_exact=FALSE,seed=12))
close_to(result$estimates,direct$estimates)
close_to(result$weights,direct$weights)
stopifnot(identical(result$data$data.IPD,sim$data.IPD),
          identical(result$data$data.AD.or,sim$data.AD.or),
          identical(design,result$simulation$design),
          identical(names(result$fits),anchored_method_names()),
          "RLB" %in% result$estimates$method)
cat("PASS: new run entry matches explicit joint-DGM generation and custom-score comparison\n")

fc <- subset(result$estimates,method=="brm-FC")
fc_bc <- subset(result$estimates,method=="brm-FC-BC" & component=="alpha")
stopifnot(all(fc$status=="available"),all(fc$interval=="Wald"),
          all(is.finite(as.matrix(fc[c("estimate","se","lower","upper","p.value")]))),
          all(is.na(fc_bc$lower)),
          identical(result$fits[["brm-FC"]],result$fits[["brm-FC-BC"]]))
close_to(fc$lower,fc$estimate+qnorm(.025)*fc$se)
close_to(fc$upper,fc$estimate+qnorm(.975)*fc$se)
cat("PASS: FC Wald reuses FC fit and returns intervals when BC is disabled\n")

scalar <- suppressWarnings(run.anchored(bc_target="coefficient",score_standardization="fixed",design,600,800,seed=12,
  score_spec=scores,effect_modifiers=character(0),return_full=FALSE))
stopifnot(identical(dim(scalar),c(5L,length(anchored_method_names()))),
          identical(colnames(scalar),anchored_method_names()),
          all(is.finite(scalar[,"brm-FC"])))
bad(get.compare.anchored(bc_target="coefficient",score_standardization="fixed",sim$data.IPD,sim$data.AD,param="RD",run_exact=FALSE))
bad(run.anchored(bc_target="coefficient",score_standardization="fixed",design,601))
bad(run.anchored(bc_target="coefficient",score_standardization="fixed",design,600,seed=NA_real_))
cat("PASS: scalar interface and invalid new-pipeline inputs\n")

batch <- suppressWarnings(simulate_anchored(bc_target="coefficient",score_standardization="fixed",design,600,800,seeds=c(12,13),
  score_spec=scores,weight_covariates="h_E",run_exact=FALSE))
first <- subset(batch$estimates,r==12,select=-r)
rownames(first) <- NULL
close_to(first,result$estimates)
summary <- summarize_anchored_simulation(batch)
stopifnot(all(is.na(summary$metrics$bias)),all(is.na(summary$metrics$coverage)),
          all(summary$metrics$n_total==2),!any(summary$metrics$component=="beta"),
          all(subset(summary$metrics,method=="brm-FC")$n_interval==2),
          all(subset(summary$metrics,method=="brm-FC-BC")$n_interval==0),
          "RLB" %in% summary$method_status$method)
## Deliberately supplied test target: verify summary arithmetic, not scientific truth.
truth <- data.frame(method="brm-FC",component="alpha",term="h_E",truth=.123)
with_truth <- summarize_anchored_simulation(batch,truth)$metrics
m <- subset(with_truth,method=="brm-FC" & component=="alpha" & term=="h_E")
d <- subset(batch$estimates,method=="brm-FC" & component=="alpha" & term=="h_E")
close_to(m$bias,mean(d$estimate)-.123)
close_to(m$coverage,mean(d$lower<=.123 & .123<=d$upper))
stopifnot(sum(is.finite(with_truth$truth))==1L)
bad(summarize_anchored_simulation(batch,transform(truth,method="removed_method")))
cat("PASS: batch retains all parameters/status, summarizes only aligned targets and requires explicit truth\n")

## BC calls remain two (BRM and FC) with EM; the extra Wald row must not refit/bootstrap.
original_exact <- exact
calls <- 0L
exact <- function(...) {
  a <- list(...); calls <<- calls+1L
  list(low=rep(-2,a$pa),up=rep(2,a$pa),p=rep(.25,a$pa), diagnostics=list(
    upper=rep(list(list(convergence=TRUE)),a$pa),
    lower=rep(list(list(convergence=TRUE)),a$pa),
    observed.null=rep(list(list(converged=TRUE)),a$pa),
    p.value=rep(list(list(diagnostics=list())),a$pa)))
}
bc <- tryCatch(suppressWarnings(run.anchored(bc_target="coefficient",score_standardization="fixed",design,600,800,seed=12,
  score_spec=scores,weight_covariates="h_E",run_exact=TRUE)),
  finally={exact <- original_exact})
stopifnot(calls==2L,identical(bc$fits[["brm-FC"]],bc$fits[["brm-FC-BC"]]))
fc <- subset(bc$estimates,method=="brm-FC" & component=="alpha")
fc_bc <- subset(bc$estimates,method=="brm-FC-BC" & component=="alpha")
close_to(fc$lower,fc$estimate+qnorm(.025)*fc$se)
stopifnot(all(fc_bc$lower == -2),all(fc_bc$upper == 2),all(fc_bc$p.value == .25))
cat("PASS: BC forwarding keeps FC Wald independent (call spy, not interval accuracy)\n")
