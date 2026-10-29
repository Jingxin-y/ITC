# RD S06 seed 9 的 BC 区间修复

只对应 `score_unit / S06 / seed 9 / brm-FC-BC`。其余逐次记录使用仓库原基础快照。
本目录 `code/` 是这次修复实际使用的完整 147 文件运行快照，逐文件匹配 `REPAIR_RUNTIME_FILES.csv`。
与基础快照只有 `RD/estimation/CI_marginal_rd_target.R` 不同；改动见 `failure_only_restart.patch`：原求解未得到有效 stationary solution 时才重启。

完整重算后的区间为 [-0.116442185313215, 0.235690918621296]，p=0.6；点估计不变。
来源说明及前后变化：[original result/repairs](../../../original%20result/repairs/RD_score_unit_S06_seed9_20261029)。
全部源文件和小型参数输入 RDS 保持原字节；逐 seed 模型 RDS 及完整拟合诊断 RDS 留在 Nibi 原压缩包中。

`rerun_bc.R` 与 `exact_recompute_arguments.rds` 是原重算入口和参数输入。显式运行以下命令会重新计算完整 BC 区间，耗时较长；本次 GitHub 同步只核对已保存结果，没有执行它。

```sh
Rscript --vanilla code/repairs/RD_stationary_restart_S06_seed9_20261029/rerun_bc.R
```
