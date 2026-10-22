## 在 anchored_v1 根目录运行；使用最新联合 RR-GOP 数据生成器。
source("data/example_anchored.R", encoding="UTF-8")
source_anchored()

weight_covariates <- c("z1","z2","z3")
model_covariates <- c("h_E","h_P")
effect_modifiers <- "h_E"
wf <- get.weight(as.matrix(data.IPD$data[,weight_covariates]),
                 data.AD$anchored.mean.x)

## 固定权重 fit_robust_lb() 已退出比较；这里仅演示保留的联合版。
## 联合 sandwich：同时计入 IPD 权重参数 gamma 和 AD 样本均值的估计误差。
rlb.joint <- fit_robust_lb_joint(data.IPD$data,wf$weight,gamma=wf$ga,
  target_mean=data.AD$anchored.mean.x,target_cov=data.AD$anchored.cov.x,
  m_AD=data.AD$anchored.n,account_for_target_mean=TRUE,
  weight_covariates=weight_covariates,model_covariates=model_covariates,
  effect_modifiers=effect_modifiers,return_full=TRUE)

rlb.joint$point.est
rlb.joint$vcov
rlb.joint$vcov.ipd
rlb.joint$vcov.target

## account_for_target_mean=FALSE 时，联合版只传播 IPD 中 gamma 的估计误差。
## account_for_gamma=FALSE 时，不再传播 gamma；同时应保持 target_mean 开关为 FALSE。
## robust 方差需要 log-binomial 成功收敛至内部解，不能解决原模型的收敛问题。
