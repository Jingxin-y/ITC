## Explicit MAIC mode preserves the weighted regression fixtures; unit mode has its own checks.
## Coefficient-BC compatibility; marginal integration is checked in check_marginal_pipeline.R.
## Legacy linear-OP regression checks; marginal defaults are tested in check_marginal_op.R.
## 在 anchored_v1 运行：Rscript estimation/tests/check_comparison_covariates.R
source("load_anchored.R", encoding="UTF-8")
source_anchored()

dat <- expand.grid(id=1:10, x=0:1, z1=0:1, z2=0:1, z3=0:1)
dat$y <- as.integer(dat$id <= 2 + dat$x + dat$z1 + dat$z3)
ipd <- list(data=dat)
ad <- get.aggre(ipd)
ad$anchored.mean.x <- c(z1=.4,z2=.6,z3=.3)
result <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", ipd,ad,run_exact=FALSE,
  weight_covariates="z2",model_covariates="z3",effect_modifiers="z1"))
stopifnot(is.list(result),identical(result$variables$weight_covariates,"z2"),
          identical(result$variables$model_covariates,"z3"),
          identical(result$variables$effect_modifiers,"z1"),
          max(abs(result$weights$weight-ifelse(dat$z2==1,1.2,.8))) < 5e-5)
for (method in c("brm","brm-BC","brm-FC-BC")) {
  fit <- result$fits[[method]]
  stopifnot(identical(unname(fit$va),unname(cbind(1,dat$z1))),
            identical(unname(fit$vb),unname(cbind(1,dat$z3))),
            identical(dim(fit$cov),c(4L,4L)),length(fit$point.est)==4L)
}

## 直接用指定的 va、vb、vg 调用底层 MLE，独立核对上层选择及其联合方差。
manual <- suppressWarnings(MLEst("IPD","RR",dat$y,dat$x,
  cbind(1,dat$z1),cbind(1,dat$z3),result$weights$weight,result$weights$ga,
  max.step=40,thres=1e-6,alpha.start=c(0,0),beta.start=c(0,0),pa=2,pb=2,
  target_mean=ad$anchored.mean.x["z2"],method="brm",
  target_cov=ad$anchored.cov.x["z2","z2",drop=FALSE],m_AD=ad$anchored.n,
  vg=matrix(dat$z2,ncol=1)))
stopifnot(max(abs(manual$point.est-result$fits$brm$point.est))<1e-10,
          max(abs(manual$cov-result$fits$brm$cov))<1e-10)
lp <- glm(y~x+z3+x:z1,data=dat,family=poisson(link="log"),weights=result$weights$weight)
stopifnot(max(abs(coef(lp)-coef(result$fits$LP)))<1e-6,
          identical(colnames(model.matrix(result$fits$LB)),c("(Intercept)","x","z3","x:z1")),
          identical(names(result$fits$RLP$point.est),names(coef(lp))),
          identical(dim(result$fits$RLP$vcov),c(4L,4L)))
cmh <- get.weighted.RR(dat,result$weights$weight,strata=dat$z3)
stopifnot(identical(cmh,result$fits$CMH))
cat("PASS: separate weight, nuisance/main-effect and effect-modifier columns propagate to all estimators\n")

tab <- result$estimates
stopifnot(sum(tab$method=="brm" & tab$component=="alpha")==2L,
          sum(tab$method=="brm" & tab$component=="beta")==2L,
          all(tab$status[tab$method %in% c("brmad","brm_b-BC")]=="not_applicable_with_effect_modifiers"),
          all(tab$component[tab$method %in% c("CMH","CMH_sandwich")]=="common_effect"),
          all(tab$interval[tab$method=="brm-FC-BC" & tab$component=="alpha"]=="BC_not_run"))

## 截距-only nuisance / 不加权 / 单个推导分数，以及默认矩阵返回的兼容性。
empty <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", ipd,ad,run_exact=FALSE,
  weight_covariates=character(0),model_covariates=character(0),effect_modifiers="z1"))
stopifnot(all(empty$weights$weight==1),ncol(empty$fits$brm$vb)==1L,
          ncol(empty$fits$brm$va)==2L,
          identical(colnames(model.matrix(empty$fits$LP)),c("(Intercept)","x","x:z1")))
dat$h_E <- (dat$z1-.3)/sqrt(.21)
score <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", list(data=dat),ad,run_exact=FALSE,
  weight_covariates="z2",model_covariates=c("z3","h_E"),effect_modifiers="h_E"))
stopifnot(identical(colnames(score$fits$brm$va),c("(Intercept)","h_E")),
          identical(colnames(score$fits$brm$vb),c("(Intercept)","z3","h_E")))
score.ad <- get.aggre(list(data=dat),covariates="h_E")
score.ad$anchored.mean.x[] <- 0
score.weight <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", list(data=dat),score.ad,run_exact=FALSE,
  weight_covariates="h_E",model_covariates="z3",effect_modifiers=c("z1","z3")))
stopifnot(identical(names(score.weight$weights$ga),"h_E"),
          ncol(score.weight$fits$brm$va)==3L,ncol(score.weight$fits$brm$vb)==2L,
          sum(score.weight$estimates$method=="brm" & score.weight$estimates$component=="alpha")==3L)
default <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", ipd,ad,run_exact=FALSE))
explicit <- suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", ipd,ad,run_exact=FALSE,
  weight_covariates=c("z1","z2","z3"),model_covariates=c("z1","z2","z3")))
stopifnot(is.matrix(default),identical(dim(default),c(5L,length(anchored_method_names()))),identical(default,explicit))
cat("PASS: full parameter output, unsupported-method labels, zero/single columns and legacy default\n")

## 真实 RR/RD 上游 BC 由 check_effect_modifiers.R / check_firth_rr_rd.R 检查。
## 这里记录调用，独立检查比较流程给 BRM/Firth 传入的矩阵、权重和参数。
## 替身输出只检查结果表装配，不是置信区间数值检验。
exact.real <- exact
bc.calls <- list()
exact <- function(...) {
  a <- list(...)
  bc.calls[[length(bc.calls)+1L]] <<- a
  list(low=rep(-1,a$pa),up=rep(1,a$pa),p=rep(.5,a$pa),diagnostics=list(
    upper=rep(list(list(convergence=TRUE)),a$pa),lower=rep(list(list(convergence=TRUE)),a$pa),
    observed.null=rep(list(list(converged=TRUE)),a$pa),
    p.value=rep(list(list(diagnostics=list())),a$pa)))
}
bc <- tryCatch(
  suppressWarnings(get.compare.anchored(model_weighting="maic",bc_target="coefficient",score_standardization="fixed",op_model="linear", ipd,ad,seed=7,
    weight_covariates="z2",model_covariates="z3",effect_modifiers="z1")),
  finally={exact <- exact.real}
)
stopifnot(length(bc.calls)==2L)
for (i in seq_along(bc.calls)) {
  fit <- bc$fits[[c("brm-BC","brm-FC-BC")[i]]]
  a <- bc.calls[[i]]
  stopifnot(identical(a$va,fit$va),identical(a$vb,fit$vb),
            identical(a$pars,fit$point.est),identical(a$se,fit$se.est),
            identical(a$weight,bc$weights$weight),
            identical(a$y,dat$y),identical(a$x,dat$x),a$pa==2L,a$pb==2L)
}
for (method in c("brm-BC","brm-FC-BC")) {
  rows <- bc$estimates[bc$estimates$method==method & bc$estimates$component=="alpha",]
  stopifnot(nrow(rows)==2L,all(rows$lower == -1),all(rows$upper == 1),
            all(rows$p.value == .5))
}
cat("PASS: selected-model BRM/Firth BC forwards full fits and assembles both alpha rows (call spy)\n")

## 标量边界修正的完整系数要与修正后的估计一致，不伪造 alpha-beta 协方差。
constant.fit <- suppressWarnings(get.estimate(op_model="linear", "AD",ipd,rep(1,nrow(dat)),model_covariates="z3"))
edge <- ipd$data
edge$y[edge$x==0] <- 0
set.seed(8)
corrected <- get.brm.ad(constant.fit,edge,rep(1,nrow(edge)))
stopifnot(corrected$coefficients[1,1]==corrected$point.est[1],
          corrected$coefficients[1,2]==corrected$se.est[1],
          corrected$cov[1,1]==corrected$se.est[1]^2,
          all(is.na(corrected$cov[1,-1])),isTRUE(corrected$bayes_adjusted))
cat("PASS: scalar Bayes full output is internally consistent and marks unavailable cross-covariance\n")
