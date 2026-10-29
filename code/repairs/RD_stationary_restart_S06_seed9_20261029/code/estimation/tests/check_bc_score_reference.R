source("load_anchored.R",encoding="UTF-8"); source_anchored()
close_to <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,scale=1)))
design <- joint_rr_gop_design(event="common",ess_ratio=.5)
scores <- make_brm_score_spec()
r <- suppressWarnings(run.anchored(design,600,800,seed=12,score_spec=scores,
  weight_covariates=paste0("z",1:6),return_data=TRUE))
set.seed(12); generated <- generate_anchored_data(600,800,design)
stopifnot(identical(r$data$data.IPD,generated$data.IPD),
          identical(r$data$data.AD.or,generated$data.AD.or),
          identical(design,r$simulation$design))
close_to(r$data$data.IPD$data$h_E,(r$data$data.IPD$data$S_E-.9)/sqrt(.63))
spec <- r$variables$score_spec
bc <- apply_brm_scores(r$data$data.AD.or,spec)$data
ac <- apply_brm_scores(r$data$data.IPD,spec)$data
close_to(unname(spec$center),c(mean(bc$S_E),mean(bc$S_P)))
close_to(unname(spec$scale),c(sd(bc$S_E),sd(bc$S_P)))
close_to(c(mean(bc$h_E),mean(bc$h_P)),c(0,0))
close_to(c(sd(bc$h_E),sd(bc$h_P)),c(1,1))
stopifnot(abs(mean(ac$h_E))>.05,abs(spec$center['h_E']-.9)>1e-5)
for (m in c("brm","brm-FC")) close_to(r$fits[[m]]$va[,"h_E"],ac$h_E)
for (m in c("LB","LP")) {
  close_to(r$fits[[m]]$model$h_E,ac$h_E)
  close_to(r$fits[[m]]$model$h_P,ac$h_P)
}
moments <- brm_score_summary(r$data$data.AD,spec)
close_to(moments$anchored.mean.x[c("h_E","h_P")],c(h_E=0,h_P=0))
close_to(moments$anchored.var.x[c("h_E","h_P")],c(h_E=1,h_P=1))
cat("PASS: true DGM scores unchanged; all analysis models share BC sample center/SD\n")

## The covariance of a sum includes correlations, not just marginal variances.
z <- data.frame(z1=c(0,0,0,1,1,1),z2=c(0,0,1,0,1,1),x=rep(0:1,3),y=rep(0:1,3))
ad <- get.aggre(list(data=z))
custom <- make_brm_score_spec(c("z1","z2"),character(0),reverse_covariates=character(0),
  reference_prob=c(z1=.5,z2=.5),data.ad=ad)
close_to(unname(custom$scale['h_E']),sd(z$z1+z$z2))
stopifnot(abs(custom$scale['h_E']-sqrt(sum(ad$anchored.var.x)))>.01)
ad$anchored.cov.x <- NULL
stopifnot(inherits(try(brm_scores_from_bc(custom,ad),silent=TRUE),"try-error"))
cat("PASS: BC score SD uses reported covariance; missing covariance is not replaced by hidden independence\n")

alpha <- design$parameters$alpha_AC + design$parameters$lambda*(spec$center['h_E']-.9)/sqrt(.63)
lambda <- design$parameters$lambda*spec$scale['h_E']/sqrt(.63)
t <- subset(r$reference_truth,method=="brm")
close_to(t$truth[match(c("(Intercept)","h_E"),t$term)],unname(c(alpha,lambda)))
bad_vars <- r$variables
bad_vars$score_spec <- brm_scores_from_bc(make_brm_score_spec(c("z1","z4")),r$data$data.AD)
stopifnot(is.null(brm_score_coefficient_truth(design,bad_vars)))
cat("PASS: per-replicate coefficient truths transform correctly; wrong EM score has no invented truth\n")

batch <- suppressWarnings(simulate_anchored(design,600,800,seeds=c(12,13),
  score_spec=scores,weight_covariates=paste0("z",1:6)))
stopifnot(length(batch$score_references)==2L,
  !identical(batch$score_references[[1]]$score_spec$center,batch$score_references[[2]]$score_spec$center))
s <- summarize_anchored_simulation(batch)$metrics
d <- subset(batch$estimates,method=="brm-FC" & component=="alpha" & term=="h_E")
t <- subset(batch$reference_truth,method=="brm-FC" & term=="h_E")
tv <- t$truth[match(d$r,t$r)]
m <- subset(s,method=="brm-FC" & term=="h_E")
stopifnot(all(d$status=="available"),m$truth_varies,m$n_coverage==2)
close_to(m$bias,mean(d$estimate-tv))
close_to(m$coverage,mean(d$lower<=tv & tv<=d$upper))
close_to(m$empirical_sd,sd(d$estimate-tv))
close_to(m$empirical_sd_estimate,sd(d$estimate))
static <- data.frame(method="brm-FC",component="alpha",term="h_E",truth=.3)
stopifnot(inherits(try(summarize_anchored_simulation(batch,static),silent=TRUE),"try-error"))
cat("PASS: simulation saves each score reference; bias/coverage/error SD align with each replicate's truth\n")
