# 非 Wald 区间共享检查入口

本目录保留当前计算源码的原字节审阅副本；实际运行仍从 `RR/code` 或 `RD/code` 加载依赖。
[SOURCE_INDEX.csv](SOURCE_INDEX.csv) 记录原文件位置及 SHA256，[METHOD_INDEX.csv](METHOD_INDEX.csv) 记录函数入口。

| 方法 | 原程序核心函数 | 审阅文件 |
|---|---|---|
| RR brm / brm-FC 的 BC（原导出标签 Exact） | `exact_marginal_target()` | [CI_marginal_target.R](RR/estimation/helpers/CI_marginal_target.R) |
| RD brm-BC / brm-FC-BC | `exact_marginal_rd_target()` | [CI_marginal_rd_target.R](RD/RD/estimation/CI_marginal_rd_target.R) |
| RD MN-ESS | `rd_weighted_arms()` → `rd_mn_interval()` → `rd_mn_score()` | [rd_methods.R](RD/RD/estimation/rd_methods.R) |

RR 的边际目标使用 log(RR)，RD 使用自然差值。BC 程序使用候选值特定的参数 bootstrap 和 Blaker 型接受区域反演；不属于百分位或 BCa 区间。
当前 coarse / confirm / refine / p-value 次数为 600 / 1000 / 1000 / 500；level=.95，thres.dicho=.001，bisection.max.step=15，正式 CI seed 为 1000000+数据 seed。
拟合权重与目标标准化权重分别传入；bootstrap 固定协变量、拟合和标准化权重，只重新生成结局。
检查端点状态、区间返回标记、边界截断和拟合失败诊断。Exact 为保留的原导出标签，不提供有限样本精确覆盖保证。

MN-ESS 使用每臂 MAIC 加权均值和 Kish 有效样本量，反演 MN score；有效计数可为小数，属于固定权重下的近似。
STC 当前结果采用 Wald 区间，具体实现和两种不确定性范围见 [STC 说明](../STC_RR/README.md)，不在本非 Wald 清单另列。

2026-10-29 修复：RD `score_unit / S06 / seed 9 / brm-FC-BC` 的实际 BC 核心为 [修复版 CI_marginal_rd_target.R](../repairs/RD_stationary_restart_S06_seed9_20261029/code/RD/estimation/CI_marginal_rd_target.R)，差异见 [failure_only_restart.patch](../repairs/RD_stationary_restart_S06_seed9_20261029/failure_only_restart.patch)。基础审阅文件仍对应其余记录。
