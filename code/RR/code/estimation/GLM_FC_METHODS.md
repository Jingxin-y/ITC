# LB / LP 的 FC 版本

2026-09-30 新增 `anchored_glm_fc.R`。主比较、`run.anchored()`、
`simulate_anchored()` 和集群全方法输出已接入下列六个方法。
统一方法名称并合并重复 BRM 后，主比较返回 18 个方法/推断标签；集群返回 18 条方法/CI 记录。

| 方法 | 点估计修正 | 系数方差 | CI |
|---|---|---|---|
| `LB-FC` | Bernoulli log-link 的均值偏倚修正；单位权重时为 Firth AS_mean | 固定解析权重的模型期望型一阶 sandwich | Wald |
| `RLB-FC` | 与 `LB-FC` 相同 | observed bread + 经验 meat；MAIC 拟合时加 gamma / AD 均值联合项 | Wald |
| `LP-FC` | Poisson log-link 的均值偏倚修正；单位权重时为标准 Poisson-Firth | 固定解析权重的 Poisson 模型型一阶 sandwich | Wald |
| `RLP-FC` | 与 `LP-FC` 相同 | 经验 sandwich；MAIC 拟合时加 gamma / AD 均值联合项 | Wald |
| `LP-BFC` | Poisson 估计方程按 Bernoulli 结局方差推导的均值偏倚修正 | 固定解析权重的 Bernoulli 模型型一阶 sandwich | Wald |
| `RLP-BFC` | 与 `LP-BFC` 相同 | 经验 sandwich；MAIC 拟合时加 gamma / AD 均值联合项 | Wald |

`BFC` 中的 B 表示 Bernoulli，不表示 BC 试验或 bootstrap CI。
六个方法均保留原 GLM 的主效应、协变量和治疗交互项。
FC 对整个回归系数向量联合求根，并不只修正治疗系数。
LB 的 log link 是非 canonical link；AS_mean 不能机械替换为最大化
`log likelihood + .5 log det(I)`。标准 Poisson log-link 是 canonical link，
单位权重时才有相应的 Jeffreys 惩罚等价性。

## 二元结局下 LP 的两种口径

`LP-FC` 以 Poisson 模型推导，保持传统 Poisson-Firth 的含义。
二元结局的真实方差是 `mu*(1-mu)`，所以 Poisson 推导不能自动保证此时的
均值偏倚消除。给它换用经验 sandwich 只改变方差，不改变点估计的修正模型。

`LP-BFC` 仍解 log-Poisson 的估计方程，但用实际 Bernoulli 方差计算
领先阶均值偏倚修正。这个版本是本项目针对二元结局的估计方程修正，
不冒称为标准 Poisson-Firth。它要求 log-risk 均值模型正确且概率在内部。

## 固定解析权重公式

记设计向量为 `X_i`，`mu_i=exp(X_i' theta)`。权重在求解时归一化为均值 1；
一个加权观察仍是一个独立结局。所有权重共同缩放不改变根或方差。

对于 LB：

```
U = sum_i w_i X_i (y_i-mu_i)/(1-mu_i)
J = sum_i w_i mu_i/(1-mu_i) X_i X_i'
B = sum_i w_i^2 mu_i/(1-mu_i) X_i X_i'
V = J^{-1} B J^{-T}
L_i = X_i' V X_i; K_i = X_i' J^{-1} X_i
A = sum_i X_i [ .5 w_i mu_i/(1-mu_i) L_i
                + w_i mu_i^2/(1-mu_i)^2 (L_i-w_i K_i) ]
```

该式是已有 BRM 固定权重 Bernoulli 偏倚修正公式在 log-link 下的化简。
单位权重时 `V=J^{-1}`，第二项消失，`A=.5 sum_i h_i X_i`。

对于 LP：

```
U = sum_i w_i X_i (y_i-mu_i)
J = sum_i w_i mu_i X_i X_i'
B = sum_i w_i^2 v_i X_i X_i'
V = J^{-1} B J^{-T}
A = .5 sum_i w_i mu_i X_i (X_i' V X_i)
```

标准 Poisson-Firth 使用 `v_i=mu_i`；Bernoulli 版本使用 `v_i=mu_i*(1-mu_i)`。
两者均解 `U+A=0`。由于 LP 分数的参数导数不依赖结局，领先偏倚中的随机
一阶导数项为零，保留二阶导数与 `V` 的收缩项。

这些修正针对回归系数的领先阶均值偏倚。取指数后的 RR、非线性边际 RR
并不自动具有相同的偏倚消除性质；估计 MAIC gamma 的额外点估计偏倚也未纳入。

## 方差和目标

普通 FC 版本始终把拟合权重视为给定，使用上式 `V`，与普通 LB / LP
的固定权重口径对应。稳健 FC 使用未修正分数的一阶 influence function，
meat 为经验分数外积。2026-10-05 起，`RLB-FC` 在 FC 点估计处使用与 `RLB`
相同的 observed bread：

```
J_obs = sum_i w_i mu_i (1-y_i)/(1-mu_i)^2 X_i X_i'
s_i = w_i X_i (y_i-mu_i)/(1-mu_i)
V_IPD = J_obs^+ [sum_i s_i s_i'] J_obs^{+T}
```

`+` 表示 Moore-Penrose 逆，与现有 `RLB` 的处理一致；不回退到 expected bread。
校准 gamma / AD 均值传播以及集群的 joint marginal Wald 也使用该 observed bread。
全 1 拟合权重时，marginal Wald 用更新后的系数协方差作 delta-method 传播。
FC 点估计和偏倚修正公式中的 `J`、`B`、`V` 仍保持上述模型期望口径。
`LB-FC` 的模型方差和 Poisson 稳健版本保持原口径；Poisson 分数的 observed
与 expected bread 相同。这些方差仍是一阶近似，不使用修正方程 Jacobian。

`model_weighting="unit"` 是主流程默认值，所有 FC 模型使用全 1 权重，
稳健版本关闭 gamma / AD 均值项。`model_weighting="maic"` 使用 MAIC 权重；
稳健版本通过校准方程传播 gamma，并加入 `Cov(Z_BC)/n_BC` 的目标均值项。
分数中心、尺度和本次设计矩阵仍条件于已构造的值。

主比较表报告 GLM 系数及 95% Wald CI。集群 `maic_wald_rows()` 另计算
MAIC 标准化的 BC 人群 A/C 边际 log RR 和 delta-method Wald CI。
unit 模式固定标准化权重；MAIC 模式的稳健 FC 额外传播权重、经验标准化
和 AD 均值。普通 FC 的边际方差仍采用固定权重 delta-method。
新增 GLM 方法不调用 BRM 的 RR+OP Exact 求解器。

## 求解、失败与使用

求解采用带回溯的 Newton 法，使用数值 Jacobian。
收敛要求最大标准化修正分数小于 `fit_control$fc_thres`，默认 `1e-6`；
最多 `fit_control$fc_max_step` 次，默认 300。有限 MLE 仅作为初值。
零事件的两臂模型可得到有限根。全成功、概率边界、奇异设计或复杂约束
可能没有可接受的内部根；未收敛结果保留到末次迭代档案，主表点估计/CI 为 NA。
Bernoulli 版本不截断概率来伪造有效根；标准 Poisson 版本可以预测大于 1，
与原 LP 的工作均值模型一致。
2026-10-07 起，取消 `LB-FC/RLB-FC`、`LP-BFC/RLP-BFC` 在边际标准化前的
逐人反事实概率硬筛选。所有 GLM 方法均保留可计算的未截断 plug-in 点估计，
即使反事实预测等于或大于 1；标准化均值仍须有限且为正。
数值可计算不等于风险概率有效。只发生方差失败时仍可恢复可计算的边际点估计。
事实拟合的 Bernoulli 概率约束、收敛判断、拟合方程、系数和方差公式保持原实现。
完整规则见 [反事实风险导出说明](COUNTERFACTUAL_RISK_CHECKS.md)。

```r
source("load_anchored.R", encoding="UTF-8")
source_anchored()
source("estimation/example_glm_fc.R", encoding="UTF-8")
result <- example_glm_fc()
result$fits[["LB-FC"]]$point.est
result$fits[["RLP-BFC"]]$vcov

# 独立拟合：dat 是含 y、x 和分析协变量的 data.frame
lb_fc <- fit_logbinomial_fc(dat, model_covariates=c("h_E","h_P"), effect_modifiers="h_E")
lp_fc <- fit_lp_fc(dat, model_covariates=c("h_E","h_P"), effect_modifiers="h_E")
lp_binary_fc <- fit_lp_fc(dat, bias_model="bernoulli",
                         model_covariates=c("h_E","h_P"), effect_modifiers="h_E")
```

集群脚本同步修正了此前 unit 点拟合、MAIC 方差/Exact 拟合混用的问题。
通过环境变量 `MAIC_MODEL_WEIGHTING=unit`（默认）或 `maic` 显式选择。
新增 FC 方法时 checkpoint settings 版本变为 4；当前版本为 10，
记录取消反事实概率硬筛选后的数值检查规则。使用新结果目录，已有历史结果不会重算。

## 核验与来源

```
Rscript estimation/tests/check_glm_fc.R
Rscript estimation/tests/check_glm_fc_pipeline.R
Rscript estimation/tests/check_glm_fc_cluster.R
Rscript estimation/tests/check_rlb_fc_observed_bread.R
Rscript estimation/tests/check_counterfactual_risks.R
```

核验包含两臂闭式解、零事件、非等权与缩放、独立数值随机 Taylor oracle、
Poisson Jeffreys 惩罚梯度、LB 单位权重 AS_mean、独立 sandwich / 联合矩阵、
系数和边际输出、真值对齐以及失败归档。小规模核验验证实现，不承诺所有
场景都减少偏倚或达到 95% 覆盖率。

标准 GLM 均值偏倚修正参考 [brglm2 作者文档](https://ikosmidis.r-universe.dev/brglm2/doc/manual.html)
和 [Kosmidis, Kenne Pagui & Sartori (2020)](https://doi.org/10.1007/s11222-019-09860-6)。
估计方程偏倚展开参考 [Kosmidis & Lunardon](https://arxiv.org/abs/2001.03786)。
固定解析权重以及 Bernoulli LP 的上述模型期望化简由本项目推导，并经独立
数值 oracle 核验；它们不是把该论文的经验偏倚修正公式直接照搬。

2026-09-30 本地 R 4.3.2 核验通过：以上三个新增测试，以及
`check_unit_model_weights.R`、`check_joint_pipeline.R`、
`check_comparison_covariates.R`、`check_failure_handling.R`、
`check_iteration_archive.R`、`check_marginal_pipeline.R`、
`check_firth_fixed_weights.R`。其中 marginal pipeline 包含实际小预算 bootstrap CI。
cluster 测试用调用记录器核验权重和导出，不把记录器输出解释为 CI 精度检验。
新增实现只依赖 base R 与项目已有 MASS。此机器未安装 `brglm2` / `sandwich`，
未声称运行了这两个包的拟合对照或依赖它们的全部旧测试。
