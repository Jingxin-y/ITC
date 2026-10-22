## A fully crossed fixture gives an independent BC-target risk-ratio oracle.
## It has effect modification and equal treatment allocation in every stratum.
source("load_anchored.R", encoding="UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R", encoding="UTF-8")

dat <- expand.grid(id=1:100, x=0:1, z1=0:1, z2=0:1)
events <- ifelse(dat$x==1, ifelse(dat$z1==1,60,20), ifelse(dat$z1==1,20,10))
dat$y <- as.integer(dat$id<=events)
ipd <- list(data=dat)
target <- expand.grid(i=1:100, j=1:100)
target$z1 <- as.integer(target$i<=20)
target$z2 <- as.integer(target$j<=60)
ad <- get.aggre(list(data=target))
oracle <- log((.8*.20+.2*.60)/(.8*.10+.2*.20))
source_oracle <- log((.5*.20+.5*.60)/(.5*.10+.5*.20))
close_to <- function(a,b,tol=1e-7)
  stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))

for (mode in c("unit","maic")) {
  z <- suppressWarnings(get.compare.anchored(ipd,ad,run_exact=FALSE,
    score_standardization="fixed", weight_covariates=c("z1","z2"),
    model_covariates=c("z1","z2"), effect_modifiers="z1",
    non_effect_modifiers="z2", model_weighting=mode))
  close_to(z$fits$CMH$est,oracle,tol=1e-6)
  close_to(z$fits$CMH_sandwich$est,oracle,tol=1e-6)
  stopifnot(abs(z$fits$CMH$est-source_oracle)>.1)
  detail <- get.weighted.RR(dat,z$weights$weight,return_details=TRUE,
                             model_covariates=c("z1","z2"))
  stopifnot(nrow(detail$part$counts)==4L,!detail$part$used_correction)
  z$data <- list(data.IPD=ipd,data.AD=ad)
  rows <- maic_wald_rows(z)
  cmh <- rows$method %in% c("CMH","CMH_sandwich")
  stopifnot(all(rows$component[cmh]=="marginal"),
    all(rows$target[cmh]=="A/C marginal log RR in BC target"),
    all(rows$fitting_weights[cmh]=="MAIC"),
    all(rows$standardization_weights[cmh]=="MAIC"),
    all(rows$status[cmh]=="available"))
  close_to(rows$estimate[cmh],rep(oracle,2),tol=1e-6)
  ## Correct log-risk regressions must standardize to the same oracle target.
  for (method in c("LB","LP","RLB","RLP")) {
    row <- rows[rows$method==method & rows$ci_type=="Wald",,drop=FALSE]
    stopifnot(nrow(row)==1L,row$status=="available")
    close_to(row$estimate,oracle,tol=1e-6)
  }
  if (mode=="unit") {
    close_to(z$fits$CMH$se,z$fits$CMH_sandwich$se)
    stopifnot(grepl("fixed MAIC weights",rows$uncertainty[rows$method=="CMH_sandwich"],fixed=TRUE))
  } else {
    stopifnot(z$fits$CMH_sandwich$var.mu>0,
      rows$uncertainty[rows$method=="CMH_sandwich"]=="CMH joint sandwich")
  }
  cat("PASS:",mode,"CMH uses four MAIC-weighted strata and the BC marginal RR oracle\n")
}
