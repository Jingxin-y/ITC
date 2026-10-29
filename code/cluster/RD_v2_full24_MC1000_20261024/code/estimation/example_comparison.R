## 在 anchored_v1 根目录运行：分别选择平衡变量、模型变量及 effect modifiers。
source("data/example_anchored.R", encoding="UTF-8")
source_anchored()

weight_covariates <- c("z1", "z2", "z3")
## 所有比较方法使用同一套分析分数；这里只改估计假设，不改变生成模型。
scores <- make_brm_score_spec(
  em_covariates = c("z1", "z2", "z3"),
  non_em_covariates = c("z4", "z5", "z6"))
model_covariates <- c("h_E", "h_P")
effect_modifiers <- "h_E"
non_effect_modifiers <- "h_P"

## 若平衡的也是分数列，先提供对应的 AD 汇总：
# weight_covariates <- "h_E"
# data.AD <- get.aggre(data.AD.or, covariates=weight_covariates)

## va=[1,h_E]；vb=[1, I(S_E=1), I(S_E=2), I(S_E=3)]（四层均出现时）。
## 非 EM h_P 在当前加权人群中被 marginalize；权重仍独立选择。
## model_covariates 继续控制 LB/LP/RLP/CMH；BRM 的新 OP 仅条件于 EM。
result <- get.compare.anchored(
  data.IPD, data.AD, run_exact=FALSE, bc_target="marginal", seed=1, score_spec=scores,
  weight_covariates=weight_covariates,
  model_covariates=model_covariates,
  effect_modifiers=effect_modifiers,
  non_effect_modifiers=non_effect_modifiers,
  op_model="marginal"
)
result$variables
result$estimates
result$marginal_estimates # 三种 BC 的 marginal logRR；开启 run_exact 后显示区间和 P 值。
result$fits$brm$coefficients
result$fits$brm$cov
result$fits$brm$op_spec          # EM / 非 EM、实际观察到的 OP 格与两臂人数
result$fits$RLP$point.est
result$fits$RLP$vcov
result$fits[["brm-FC"]]$coefficients # FC 点估计、按模型权重设置计算的 SE 与 Wald 区间
result$fits$RLB$vcov     # 稳健 log-binomial；MAIC 拟合时传播 gamma 和 AD 均值误差

## 将上面 run_exact 改为 TRUE，即显式调用 exact_marginal_target()。
## result$marginal_estimates 显示 marginal logRR、区间/P 值及状态；SE 留空。
## result$bc$brm.bc / result$bc$firth.bc 保存完整结果，effect.* 为 RR 尺度。
## 旧逐 alpha 区间必须显式设置 bc_target="coefficient"。
## 含修饰效应时，CMH 仍返回合并 RR；brmad 和 brm_b-BC 标为不适用。
## Wald 行仍为条件模型参数及 CMH 合并效应；marginal BC 行为标准化的 A/C logRR。
## BRM 的条件参数现在只条件于 h_E；h_P 已被 marginalize。
## 对照旧线性 OP 时显式传 op_model="linear"，此时 vb=[1,model_covariates]。

## 一次生成并比较的等价入口：n_AC/n_BC 都是整项试验人数。
# result <- run.anchored(design=design,n_AC=1000,n_BC=1000,seed=1,
#   score_spec=scores,weight_covariates=weight_covariates,run_exact=FALSE)

## 可选：不平衡协变量、nuisance 仅截距；没有修饰变量时也可要求完整结果。
# result <- get.compare.anchored(data.IPD,data.AD,run_exact=FALSE,
#   weight_covariates=character(0),model_covariates=character(0),
#   effect_modifiers=NULL,return_full=TRUE)
