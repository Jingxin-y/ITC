# 固定 MAIC 权重下的 FC 实现（2026-09-20）

当前 RR/RD 的 `method="firth"` 和比较入口 `brm-FC-BC` 已使用固定代表性权重下的均值偏倚修正。接口名称沿用 Firth；严格说，非等权时这是对加权估计方程的模型期望型偏倚修正。它不把 MAIC 权重解释为独立重复观测的频数。

## 适用范围

给定协变量、治疗分组和非负权重，结局是相互独立且均值模型正确的 Bernoulli 观测。假设参数维度固定、信息矩阵可逆、概率处于内部、存在正则一致根且少数观测的权重不会渐近主导样本。

修正针对当前 alpha/beta 参数的一阶均值偏倚：RR 为 log RR/log OP，RD 为 atanh RD/log OP。非线性变换后的效应不自动具有相同的均值偏倚性质。这不是任意有限样本严格无偏或所有稀疏数据必有有限根的保证。

## 校正公式及推导

令 \(\theta=(\alpha^\top,\beta^\top)^\top\)，并定义

\[
 e_i=Y_i-p_i,\quad v_i=p_i(1-p_i),\quad
 D_i=\nabla p_i,\quad H_i=\nabla^2 p_i,\quad q_i=D_i/v_i.
\]

原加权分数和相关矩阵为

\[
 U=\sum_iw_iq_ie_i,\qquad
 I=\sum_iw_iD_iD_i^\top/v_i,\qquad
 B=\sum_iw_i^2D_iD_i^\top/v_i,\qquad
 V=I^{-1}BI^{-\top}.
\]

这里 \(I=-E(U')\)，\(B=\operatorname{Var}(U)\)。权重放大同一个人的随机贡献，因此分数方差使用 \(w_i^2\)，不是将该人变成 \(w_i\) 个独立观测。\(V\) 是固定权重的模型期望型一阶 sandwich。

令

\[
 C_i=\nabla q_i=H_i/v_i-(1-2p_i)D_iD_i^\top/v_i^2.
\]

则 \(U'=-I+R\)，其中 \(R=\sum_iw_ie_iC_i\)。在真参数处展开 \(U(\hat\theta)=0\)，令 \(\delta_1=I^{-1}U\)：

\[
 b_{\rm leading}
 =I^{-1}\{E(RI^{-1}U)+\tfrac12E(U''):V\}.
\]

不同人的结局相互独立，所以

\[
 E(RI^{-1}U)=\sum_iw_i^2C_iI^{-1}D_i,
\]

\[
 \tfrac12E(U''):V
 =-\sum_iw_iC_iVD_i
  -\tfrac12\sum_iw_iD_i\,\operatorname{tr}(VH_i)/v_i.
\]

在 \(U+A=0\) 中，增加 \(A\) 造成的一阶参数位移为 \(+I^{-1}A\)。因此取 \(A=-I b_{\rm leading}\)，得到代码采用的公式：

\[
 \boxed{
 A=\tfrac12\sum_iw_iD_i\,\operatorname{tr}(VH_i)/v_i
   +\sum_iw_iC_i(V-w_iI^{-1})D_i.}
\]

第二项一般非零；只把原 trace 项中的 \(I^{-1}\) 换成 \(V\) 不完整。实现通过 alpha/alpha、alpha/beta、beta/beta 三个矩阵块进行收缩，无需逐人三维 Hessian 数组。

随机 Taylor 展开和估计方程偏倚消除框架参考 [Kosmidis & Lunardon, Empirical bias-reducing adjustments to estimating functions](https://arxiv.org/abs/2001.03786)，第 3.2 节及补充材料 S3。上面的固定权重 Bernoulli 模型期望化简是本次针对代码的推导，不将其冒充文中的经验修正公式或原始 Firth 论文的非等权结论。

- 当全部 \(w_i=1\)，\(B=I\)、\(V=I^{-1}\)，第二项为零，退化为原普通 Firth 校正。
- 将所有权重乘以 \(c>0\)，\(U\) 和 \(A\) 同时乘以 \(c\)，\(V\) 不变，因此根不变。
- 一个权重为 2 的人仍只有一个 Bernoulli 结局；复制成两个独立观测是不同的抽样机制。旧整数频数复制测试已替换为上述缩放不变性测试。
- 计算 FC 时使用的是模型矩阵 \(B\)，不是把经验残差 meat 不加区分地代入偏倚公式。

## 拟合与诊断

2026-10-01：两臂 RR/OP 概率映射 `getProbRR()` 与 `getProbScalarRR()` 已共享数值稳定的精确公式，移除旧的 ±12 参数边界近似。令 `a=abs(logRR)`、`r=exp(-a)`，较大风险为 `2/(1+r+sqrt((1-r)^2+4*exp(-a-logOP)))`，较小风险为其乘以 `r`；实际计算在对数尺度完成，并用 `expm1()`、`log1p()` 避免相消和溢出。仅保留浮点舍入的概率上界保护，不用固定 epsilon 把概率强行截到内部。极端参数仍可能因机器精度返回 0/1，FC 的内部概率检查与原收敛容差保持有效。

此次修改作用于共享概率映射，因此普通 BRM 的点估计及 RR 方差也使用新映射；三臂 GOP 数据生成与 GLM-FC 方程保持原实现。集群检查点 settings.version 更新至 5，避免将使用旧概率映射生成的结果直接复用。历史结果及冻结代码快照保留原计算口径，不能视为本次修复后重新运行的结果。

`helpers/firth_augmentation.R` 实现 score、两部分 adjustment、I、B、V；`helpers/MLE_Point_Firth.R` 对所有 alpha/beta 联合求根。

求解器将权重归一化为均值 1。令初值的信息矩阵为 \(I_0\)，固定尺度 \(s_j=\max\{1,\sqrt{(I_0)_{jj}}\}\)。成功要求

\[
 \max_j |(U_j+A_j)/s_j|\leq \texttt{thres}.
\]

`get.estimate()` 默认 `thres=1e-6`，FC 最多允许 300 次 Newton 迭代，可通过 `fit_control$fc_max_step` 调整。Jacobian 使用中央数值差分，Newton 步限制长度并回溯；差分或候选参数无效时采用单侧差分或拒绝候选。概率必须处于未截断的 \((10^{-8},1-10^{-8})\) 区间。

初始信息矩阵不可逆时报错；迭代耗尽或残差不能继续下降时，底层返回 `convergence=FALSE`。上层 `MLEst()` 返回标记 `converged=FALSE` 的末次拟合，保留 `last_optimizer` 和 `optimization_input`，SE/区间为 NA。比较入口将它单独归档到 `nonconverged_fits`，不计入主汇总，也不继续计算 BC。

返回对象新增：

- `fit$fc.diagnostics`：求解方式、迭代次数、分数评估次数、退出原因、最大标准化残差、原分数残差和尺度。
- `fit$fc.weight.interpretation == "fixed_analytic"`。
- `fit$variance.method`：明确方差分支及其一阶性质。

## 方差与 BC

- `group="AD", method="firth"`：使用当前 FC 点的 \(V=I^{-1}BI^{-\top}\)，直接接受固定权重。单位权重时等于普通逆信息矩阵；非等权时不再误用逆 Fisher。
- `group="IPD", method="firth"`：沿用原联合 sandwich，在 FC 点计算，包含原流程的 gamma/AD 均值传播。偏倚修正不改变一阶渐近分布；这里不声称推导了 FC 的高阶或有限样本精确方差。原 AD 不确定性计算失败后的 warning/回退行为保留。
- IPD 分支必须传与实际权重匹配的 `gamma.ipd`/`gamma.weight`。现在会比较重建权重与点估计权重的相对值；不匹配时报错。应传递 `get.weight()` 返回的 `ga`；仅希望固定权重方差时用上述 AD 分支。
- BC 保持用户确定的固定协变量、分组及权重机制，默认调用 `exact_marginal_target()`，对 FC 参数给出的 marginal logRR 反演区间；内部使用普通加权 likelihood，optim.maxit=200，协方差仅确定搜索宽度。显式 `bc_target="coefficient"` 保留旧 `exact()`，其 max.step=100、optim.maxit=300。
- 普通 `method="brm"` 的默认迭代上限和未收敛归档规则见 `FAILURE_HANDLING_AND_SCORES.md`。

## 调用示例

```r
wf <- get.weight(covariate_matrix(ipd$data), ad$anchored.mean.x)

fit <- get.estimate(
  "IPD", ipd, wf$weight, gamma.ipd = wf$ga,
  target_mean = ad$anchored.mean.x,
  target_var = ad$anchored.var.x, m_AD = ad$anchored.n,
  target_cov = ad$anchored.cov.x,
  method = "firth", param = "RR"
)
fit$fc.diagnostics

# 任意固定代表性权重的 FC 与固定权重一阶方差
fit.fixed <- get.estimate("AD", ipd, wf$weight, method = "firth", param = "RD")
```

选择权重/模型变量或 effect modifiers 的现有参数继续适用。比较入口中的 `brm-FC-BC` 自动使用新实现，无需改方法名。

## 已运行的验证

`tests/check_firth_fixed_weights.R` 包含独立分数数值微分的随机 Taylor oracle（RR/RD 各 pa=1/2/3）、两臂加权风险率的解析 delta-method 偏倚抵消、单位权重退化、权重缩放、完整拟合的根残差、固定权重及联合方差、gamma 不匹配、零权重、迭代失败和奇异设计检查。

本次六组 Taylor 检查的最大绝对差约 \(3.9\times10^{-7}\)，两臂一阶偏倚抵消误差小于 \(10^{-10}\)。这验证推导和代码一致，不代替大规模偏倚/覆盖率实验。

额外生成 40 个 n=200 的 RR/RD 样本（各含 2 个 alpha、3 个 beta；log OP 截距为 -3 或 -6），38 个达到根容差。2 个 log OP=-6、总事件数均为 11 的样本在不同初值下均未达到容差，保留失败诊断，不声称不存在其他根。这是小规模求解器实测，不将 38/40 解释为一般收敛率或覆盖率；没有为强行得到有限结果而修改校正方程。

`check_firth_rr_rd.R`、`check_effect_modifiers.R`、`check_weight_covariates.R`、`check_comparison_covariates.R` 也已通过，包含真实 RR/RD 小规模 BC 调用。额外的 `check_robust_logbinomial.R` 在此机器因未安装 `sandwich` 包而无法执行；该脚本测试的 RLB 实现不属于本次 FC 修改。
