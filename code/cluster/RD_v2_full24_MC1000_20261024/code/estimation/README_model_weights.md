# 模型拟合权重与边际标准化权重

从 2026-09-29 起，`get.compare.anchored()`、`run.anchored()` 和
`simulate_anchored()` 默认使用 `model_weighting="unit"`：

1. 先按 `weight_covariates` 和 BC 汇总均值计算原 MAIC 权重及诊断。
2. 回归结局模型拟合使用 `rep(1, n_AC)`。
3. 现有 marginal logRR 输出仍使用原 MAIC 权重作标准化，目标保持为 BC 人群。
4. 从 2026-10-06 起，CMH / CMH_sandwich 的分层计数始终使用 MAIC 权重，
   不受回归模型的 `model_weighting` 控制。

三个 EM、三个非 EM 的数据生成机制及分析分数保持原设定。
BC 样本均值/协方差仍用于分析分数的中心与尺度；这一步不是模型加权。

## 计算方式

结局模型的似然/估计方程对每个 IPD 观察取相同权重。
得到拟合风险 `p_A_hat(z_i)`、`p_C_hat(z_i)` 后，已有的边际效果仍为：

```text
log_RR_target = log(sum(w_MAIC * p_A_hat) / sum(w_MAIC * p_C_hat))
```

BRM、FC、LB、LP、稳健 LB/LP，以及可用时的标量 Bayes 修正，
均使用全 1 的模型权重。BC 内部的受约束似然和模拟后的模型拟合也使用全 1 权重。
2026-09-30 新增的 LB / LP FC 方法同样遵循此设置，详见 `GLM_FC_METHODS.md`。
边际 BC 的 `standardization.weight` 单独传入 MAIC 权重。

模型系数的方差不再传播 MAIC gamma 或校准目标均值的估计误差：
BRM/FC 使用已有的固定权重方差路径，稳健 LB/LP 与 CMH 关闭权重估计的联合方差项。
现有边际 BC 仍条件于给定的 MAIC 标准化权重，未新增重估 MAIC 的 bootstrap。
BC 分数标准化仍沿用既有处理方式。

CMH 仍按 `model_covariates` 的组合分层，并用 MAIC 加权计数合并。
在当前随机试验设计的各层治疗分配概率相同、且校准足以恢复目标分布时，
其大样本目标为 BC 人群边际 RR；有限样本的合并值不必等于回归标准化值。
`unit` 模式下，CMH 和 CMH_sandwich 都固定 MAIC 权重，使用原分层计数影响函数方差，
两行数值相同；`maic` 模式下，只有 CMH_sandwich 额外传播 gamma 和 BC 目标均值。
普通 CMH 始终使用固定权重方差。该开关保持与其他稳健方法的既有不确定性范围一致。
CMH 固定权重方差对经验分层计数求影响函数，而回归 delta 方差固定标准化协变量；
因此不能把两者称为完全相同的条件方差。

## 接口与记录

```r
source("load_anchored.R", encoding="UTF-8")
source_anchored()
design <- webster_rr_design(event="common", ess_ratio=.5)
result <- run.anchored(design, n_AC=1000, n_BC=1000, seed=1,
  score_spec=make_brm_score_spec(), weight_covariates=paste0("z",1:6),
  model_weighting="unit", run_exact=FALSE)

result$weights$weight       # 已计算的 MAIC 权重，用于边际标准化
result$model_weights$weight # 回归模型权重，unit 模式全部为 1
result$variables$model_weighting
result$variables$cmh_weighting # 始终为 maic
result$variables$cmh_sandwich_uncertainty # unit 为 fixed；maic 为 joint
result$marginal_estimates
```

矩阵返回模式分别用 `attr(result,"weights")` 和 `attr(result,"model_weights")`。
批量结果记录 `model_weighting`；逐次诊断中 `ess` 仍为 MAIC ESS，
新增 `model_ess` 是实际拟合权重的 ESS（unit 模式等于 n_AC）。
未收敛模型的记录也同时保存两套权重。

MAIC 校准失败时，未加权模型仍可拟合，条件系数按其真实拟合状态报告；
需要 MAIC 权重的边际效果标为 `weighting_failed`，不伪装成 IPD 人群效果。

若要显式复现以前的加权模型，三个入口均支持 `model_weighting="maic"`。
主配置已明确使用 `"unit"`。通用底层函数如 `get.estimate()`、`fit_lp()`
继续尊重调用者显式传入的权重，便于独立验证和旧实验复现。
独立调用 `get.brm.bc()` 也默认 `model_weighting="unit"`，传入的
`brm.est` / `bayes.est` 应来自全 1 权重拟合。原 `weight=w_MAIC` 默认作为
目标标准化权重保留；也可显式传 `standardization.weight=w_MAIC`。
候选 psi 的 nuisance likelihood 及 bootstrap likelihood 使用全 1，
而 marginal 风险平均与每次 RR 截距求根均使用 MAIC 目标权重。
旧的加权调用需显式传 `model_weighting="maic"`，并提供相应加权模型拟合。
`example_marginal_exact()` 对 BRM 和 FC 都使用同一套全 1 拟合权重，
通过共享 adapter 传递 MAIC 目标权重，`ci_control` 不能覆盖两套权重。

结果表的统计字段保持原格式，方法标签采用统一的新名称；`brm` / `RLP` 等名称不表示
实际采用 MAIC 拟合，应结合 `model_weighting` 解释。
已有结果不自动重算。建议新结果写入单独目录，例如 `results_unit_fit_maic_target`。

集群输出中 CMH 的 `component="marginal"`，填写同一 BC marginal population truth，
`fitting_weights` / `standardization_weights` 均记录为 `MAIC`；底层参数表保留
`common_effect` 标签表示合并统计量，不表示 EM 斜率。checkpoint settings 当前为版本 10，
拒绝复用旧版本结果；已有结果不自动重算，需使用新的结果目录。

2026-10-07 起，取消 Bernoulli FC 的标准化前反事实概率硬筛选，保留可计算的
未截断边际点估计。unit 与 maic 均采用这一规则，标准化均值须有限且为正；
详见 [反事实风险导出说明](COUNTERFACTUAL_RISK_CHECKS.md)。

检查：`Rscript estimation/tests/check_cmh_maic_target.R`、`check_unit_model_weights.R`；
`check_marginal_pipeline.R` 包含小规模实际 BC 运算，仅验证实现，不验证覆盖率。

## 方法名称（2026-09-30）

主比较只保留一个 `brm`，使用所选的 `model_weighting`。
所有方法名去掉 `_an`，稳健 log-binomial 统一为 `RLB`。
anchored 由代码和输出文件名中的 `anchored` 标识；结果元数据另记录 `comparison_type="anchored"`。
方法名称与迁移对应见 [METHOD_NAMES_ANCHORED.md](METHOD_NAMES_ANCHORED.md)。
