## 从 anchored_v1 根目录运行；只自定义估计分数，复用同一份真实生成数据。
source("data/example_anchored.R",encoding="UTF-8")
source_anchored()

## 只需修改这两组原始协变量；默认反向变量为选中的 z1、z3、z5。
scores <- make_brm_score_spec(
  em_covariates=c("z1","z4"),
  non_em_covariates=c("z2","z3","z5","z6")
)
## 当前示例故意选错 EM。正确组成可设为 make_brm_score_spec()。
## 如需全部正向相加，在构造器中增加 reverse_covariates=character(0)。
## 构造器先保存变量选择；比较入口默认用 BC 提供的样本矩重算中心和尺度。
scores

result <- get.compare.anchored(
  data.IPD,data.AD,score_spec=scores,
  weight_covariates=paste0("z",1:6),run_exact=FALSE
)
result$fits$brm$coefficients
result$fits$brm$op_spec$cells
result$variables$score_spec

## 单独 BRM/Firth 也可传同一配置；权重变量仍独立指定。
scores <- brm_scores_from_bc(scores,data.AD)
weight_cols <- paste0("z",1:6)
wf <- get.weight(covariate_matrix(data.IPD$data,weight_cols),
                 data.AD$anchored.mean.x[weight_cols])
fit <- get.estimate(
  "IPD",data.IPD,wf$weight,wf$ga,method="firth",score_spec=scores,
  weight_covariates=weight_cols,target_mean=data.AD$anchored.mean.x,
  target_cov=data.AD$anchored.cov.x,m_AD=data.AD$anchored.n
)
fit$coefficients

## 若只按自定义 h_E 求权重，比较入口会从原始 AD 汇总转换目标矩：
# result_score_weight <- get.compare.anchored(
#   data.IPD,data.AD,score_spec=scores,weight_covariates="h_E",run_exact=FALSE)
## 手工求分数权重时：
# ipd_est <- apply_brm_scores(data.IPD,scores)
# ad_est <- brm_score_summary(data.AD,scores,"h_E")
# wf_score <- get.weight(covariate_matrix(ipd_est$data,"h_E"),
#                        ad_est$anchored.mean.x["h_E"])
## 将 ipd_est、wf_score 及 ad_est 中的目标矩一起传给 get.estimate()。
## 完全 marginal 的 BRM：显式 effect_modifiers=character(0)。
## 默认 OP 中 h_P 被积分掉；仅改 h_P、保持 h_E 和权重不变，不会改变 BRM 点估计。
