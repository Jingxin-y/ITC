# 非 Wald 置信区间：共享检查入口

这里单独列出当前结果中的 BC 和 MN-ESS 区间。`RR/`、`RD/` 中是主程序源码的原样审阅副本；文件哈希和函数起始行见 `METHOD_INDEX.csv`、`SOURCE_INDEX.csv`。运行时仍从包根目录的 `RR/code` 或 `RD/code` 加载完整依赖。

## 方法清单

| 尺度 | 对应方法 | 区间 | 主要函数 | 审阅源码 |
|---|---|---|---|---|
| RR | brm | BC，导出标签为 `Exact` | `exact_marginal_target()` | `RR/estimation/helpers/CI_marginal_target.R` |
| RR | brm-FC | BC，导出标签为 `Exact` | `exact_marginal_target()` | 同上，传入 FC 拟合的目标点估计 |
| RD | brm（导出方法名 `brm-BC`） | BC | `exact_marginal_rd_target()` | `RD/RD/estimation/CI_marginal_rd_target.R` |
| RD | brm-FC（导出方法名 `brm-FC-BC`） | BC | `exact_marginal_rd_target()` | 同上，传入 FC 拟合的目标点估计 |
| RD | MN_eff | MN-ESS | `rd_weighted_arms()` → `rd_mn_interval()` → `rd_mn_score()` | `RD/RD/estimation/rd_methods.R` |

## RR / RD 的 BC

目标是 BC 人群中标准化的边际 A/C 效应。此处“BC 人群”指试验 BC 的目标总体；区间标签 BC 指下面的计算程序。

- RR：`psi = log(sum(v*p_A) / sum(v*p_C))`，求解器的 estimate/low/up 在 log(RR) 尺度；`effect.*` 字段为指数变换后的 RR。
- RD：`psi = sum(v*(p_A-p_C))/sum(v)`，estimate/low/up 均在自然 RD 尺度。
- `weight` 是模型拟合权重；`standardization.weight` 是目标人群标准化权重。两者应分别核对。

程序使用候选值特定的参数 bootstrap 和 Blaker 型 acceptability 反演：

1. 固定观察数据的目标点估计为 anchor；对每个候选 psi 拟合约束模型。
2. 从该候选模型的 Bernoulli 风险生成 bootstrap 结局。
3. 每个 bootstrap 样本分别拟合“候选 psi”与“固定的观察目标 anchor”，计算二者的拟合准则差。
4. 用双尾排序构造 acceptability；通过向外搜索、确认和二分细化得到 95% 接受区域的上下端点。

它不是百分位 bootstrap 或 BCa，也不等同于每次都重新估计 unrestricted MLE 的常规 bootstrap LR。FC 改变输入拟合与观察目标，随后使用同一套区间反演逻辑。协变量、治疗分组、拟合权重和标准化权重固定，仅重新生成结局；没有在 bootstrap 中重新估计 MAIC 权重或重新抽样 AD 均值。

### 正式参数

| 参数 | RR / RD 默认值 |
|---|---|
| level | 0.95 |
| coarse / confirm / refine / p-value bootstrap 次数 | 100 / 200 / 200 / 500 |
| thres.dicho | 0.001 |
| bisection.max.step | 15 |
| CI seed | 正式流程为 1000000 + simulation seed |
| 数值优化 | 原尝试优先；失败后多初值与 nlminb 后备 |

对应结果中个别 RR 修复调用把二分上限显式设为 16；源码全局默认仍为 15。若核对某次已保存区间，应查看该对象的实际 settings，不能只看函数默认值。`psi.se` / `vcov.pars` 用于确定搜索宽度，不是 bootstrap 区间提供的标准误。

检查 `interval.available`、`status`、`p.status`、端点状态、`failed.profile.fits`、`failed.bootstrap.draws` 和搜索范围。`boundary_truncated` 表示搜索边界截断；失败 bootstrap 不会被静默丢弃。`Exact` 是保留的导出标签，plug-in nuisance 参数与有限 bootstrap 次数不提供有限样本精确覆盖保证。

### 主程序调用链

- RR：`maic_run_seed()` / `marginal_exact_from_fit()` → `exact_marginal_target()`；适配层见 `RR/estimation/helpers/marginal_comparison.R`。完整批量入口在包根目录 `RR/code/estimation/cluster_all_methods.R`。
- RD：`get.compare.anchored.rd()` → `rd_marginal_bc_from_fit()` → `exact_marginal_rd_target()`；适配层见 `RD/RD/estimation/rd_comparison.R`。

## RD 的 MN-ESS

分别用 MAIC 加权的两臂结局均值估计 p_A、p_C，并在每臂计算 Kish 有效样本量 `n_eff=(sum(w))^2/sum(w^2)`。针对候选差值 d，约束 `q_A=q_C+d`，估计 nuisance 风险 q_C，再用带有限样本修正的 MN score 统计量与 `qchisq(.95,1)` 的交点构造区间。根搜索范围是 [-1,1]，数值精度 `1e-9`。

这是使用固定 MAIC 权重与可能为小数的有效计数所作的 MN 近似，不是整数二项样本的精确区间。它不提供模型 SE；`se_from_ci_width` 只是区间宽度诊断，绘图的 Accuracy 不把它当作模型标准误。

## 运行依赖

先执行 `source_anchored()` 或 `source_anchored_rd()` 加载完整代码，再通过上述适配层传入已经拟合的对象。单独 source 审阅文件不代表依赖齐全。检查源码与加载不需要运行全部 bootstrap；本次打包没有重算完整 BC。
