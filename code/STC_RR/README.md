# RR STC：直接复用已有拟合结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Date: 2026-10-26
- Version: `STC_RR_reuse_20261026_v1`
- 状态：独立后处理代码；使用保存的 RDS 完成局部验证，见 `VALIDATION.md`。未运行 Nibi 的全部 8,000 份生产 RDS。

## 做了什么

这个程序读取原来 `score_unit` 的逐 seed RDS，复用已经拟合好的模型系数、完整协方差、分数编码和模型设计，只把原来的 MAIC 目标平均换成 STC 目标平均。

**不会重新生成数据、重新拟合模型、重新做 FC 校正或计算模拟反演区间。原代码、原 RDS 和原 MAIC 结果均只读。**

运行仅依赖标准 R（`base`、`stats`、`utils`、`tools`），不需要原 simulation 代码目录、`MASS`、`sandwich` 或 `brmplus`。必要的 BRM 风险公式和分数变换已放在独立的 `stc_reuse.R` 内；没有加载生成器和拟合器。

## 输入

默认原结果根目录：

```text
/scratch/yanjin41/RR_v2_full24_MC1000_results_20261024
```

读取八个情形的：

```text
cells/S01_common_n100_ess015_score_unit/by_seed/seed_000001.rds
...
cells/S08_rare_n500_ess050_score_unit/by_seed/seed_001000.rds
```

| 情形 | event | 每项试验总人数 | 名义 Full ESS 比例 |
|---|---|---:|---:|
| 1 | common | 100 | 0.15 |
| 2 | common | 100 | 0.50 |
| 3 | common | 500 | 0.15 |
| 4 | common | 500 | 0.50 |
| 5 | rare | 100 | 0.15 |
| 6 | rare | 100 | 0.50 |
| 7 | rare | 500 | 0.15 |
| 8 | rare | 500 | 0.50 |

只接受原始 `scores + unit`、`diffuse / high / all`、`webster2024_rr_ess_em_v2`、风险目标版本 `BC_CAB_20261015`。逐文件检查 cell、seed、样本量和保存的设置是否一致。

所需 RDS 内容：

- `result$fits`：BRM/BRM-FC 的 `point.est/cov/va/vb`，普通 GLM 的系数、terms、QR，稳健/FC wrapper 保存的 `vcov`。
- `result$variables$score_spec`：本次拟合实际使用的分数定义。
- `result$data$data.AD`：BC 汇总均值、协方差和人数。
- `result$data$data.IPD$data` 中的 AC 协变量：只用于核对并映射 BRM 原设计行，不使用其结局重新估计。
- `summary`：读取原 MAIC Wald 结果用于配对比较，并且仅从这里读取原目标真值用于评价误差。

本地下载的 `RR_score_unit_combined_results_20261026.tar.gz` 只有估计结果表，**未包含上述完整拟合对象**。系数 CSV 也不提供完整协方差、分数和 nuisance 格信息，不能单独完成本任务。应在 Nibi 原 RDS 所在位置直接运行，或另行复制这些 RDS 后在本地运行。

## Nibi 运行

将新代码包上传、解压到独立目录。下面只上传本次新包：

```powershell
$pkg = "D:\codexresults\ITC\rorqual\compareITC\anchored\STC_RR_reuse_20261026.tar.gz"
scp "$pkg" "${pkg}.sha256" yanjin41@nibi.alliancecan.ca:~/
```

在 Nibi 上：

```bash
cd "$HOME"
sha256sum -c STC_RR_reuse_20261026.tar.gz.sha256
tar -xzf STC_RR_reuse_20261026.tar.gz
cd STC_RR_reuse_20261026
sha256sum -c SHA256SUMS
module load StdEnv/2023 gcc/12.3 r/4.3.1

Rscript --vanilla run_stc_reuse.R \
  --input-root /scratch/yanjin41/RR_v2_full24_MC1000_results_20261024 \
  --output-root /scratch/yanjin41/STC_RR_reuse_results_20261026
```

这条命令处理 8 个情形 × 1,000 个已有 seed，共 8,000 个来源位置。每个位置固定输出 12 个模型 × 2 种区间范围 = 24 行，共 192,000 行；失败也保留对应行。这里的 1,000 是读取原来的 1,000 个 seed，不是新做 1,000 次 simulation。

可先选一个已有 seed 检查生产 RDS 接口，使用另一个输出目录：

```bash
Rscript --vanilla run_stc_reuse.R \
  --input-root /scratch/yanjin41/RR_v2_full24_MC1000_results_20261024 \
  --output-root /scratch/yanjin41/STC_RR_reuse_smoke_20261026 \
  --scenarios 2 --seeds 17
```

`--scenarios` 与 `--seeds` 支持 `1:8`、`1,3,5`、`1:10,17`。`--level` 默认 0.95。抽取子集时汇总分母相应为实际请求的 seed 数；生产默认每情形分母仍是 1,000。

程序逐 seed 保存派生 RDS。相同命令再次运行时，只有源文件 MD5、STC 代码 MD5 和配置完全一致才复用派生结果。缺失源文件不会缓存成永久成功；补齐原 RDS 后再次运行会重新读取。不同子集、版本或配置请换输出目录。

输入、输出根目录不能相同，也不能互相包含；写文件前还会解析已有的符号链接/junction。输出使用独立配置和运行锁。若作业被强制终止后遗留 `.stc_reuse_lock`，先确认没有同目录作业运行，再用 `rmdir /scratch/yanjin41/STC_RR_reuse_results_20261026/.stc_reuse_lock` 移除这个空锁目录。

退出码：0 表示来源记录均通过入口检查；2 表示存在缺失、损坏、未完成或不兼容的来源 RDS，报告已保存。单个模型原先拟合失败不会触发重拟合；仍需查看方法返回率，不能把退出码 0 当作所有模型区间均成功。

## 方法和不确定性范围

目标保持为 **BC 总体中的边际 A/C log RR**：

```text
psi = log(sum_z pi_BC(z) * p_A_hat(z))
    - log(sum_z pi_BC(z) * p_C_hat(z))
```

根据每次保存的 BC 样本均值重构六个独立 Bernoulli 的联合分布，精确枚举 64 个二元组合；没有伪个体抽样误差。**独立性是明确的工作假设，在当前 DGP 下成立；单凭均值并不能识别一般的联合分布。** 不读取真实总体概率、风险表或 BC 个体协变量用于估计。

模型：`brm`、`brm-FC`、`LB`、`LP`、`RLB`、`RLP`、`LB-FC`、`RLB-FC`、`LP-FC`、`RLP-FC`、`LP-BFC`、`RLP-BFC`。这里 LP/RLP 指原代码的 log-Poisson/稳健 log-Poisson。CMH 保留为原 MAIC 方法，不转换为 STC；当前不适用的 `brmad` 不加入。

- `fixed_bc_summary`：固定保存的 BC 汇总量和分数编码，使用完整系数协方差作 delta 推断。
- `bc_mean_delta`：在上项基础上增加 `g_mu' * (Sigma_BC / n_BC) * g_mu`。保留 BC 报告协方差的全部非对角项。

两种区间使用相同点估计。后者只额外传播**用于目标分布重构的 BC 均值**的不确定性；仍固定模型结构和本次分数编码，不是全流程推断，不重新估计 score 的中心/尺度或协方差，也不评估分布假设错误。FC 沿用原 RDS 保存的协方差及其一阶近似。

对于 `LB/LP`，`stats::vcov()` 从已有 QR 和离散参数恢复原常规协方差，不拟合模型。稳健和 FC 版本直接使用各自保存的 wrapper `vcov`；缺失时不会用普通 GLM 协方差替代。

BRM 保留原 reduced model：nuisance 在观察到的 EM 分数格中饱和，非 EM 已 marginalize。利用原始 `S_E` 映射并复用保存的 `va/vb`，不重排 nuisance 参数，也不补格。其 `risk_C/risk_A` 是 AC reduced-model 风险在 BC EM 分布下的平均，不能直接解释成 BC 绝对风险。当前独立、乘法风险 DGP 下非 EM 共同因子在 RR 中抵消，目标 RR 保持不变。结果表包含 `risk_interpretation`。

## 失败处理

- 原拟合未收敛：`fit_not_converged`，不修复、不重拟合。
- 正概率目标 EM 格没有 nuisance 参数：`target_support_missing`，不删除该格后归一化。
- 系数协方差缺失或不合法：保留可计算的点，区间 `se_failed`。
- BC 协方差缺失：固定汇总量区间仍可用，均值 delta 区间单独失败。
- BC 均值在 0/1 边界：保留退化分布的点/可计算的固定区间，均值 Wald 标记 `target_mean_boundary`，不加伪计数。
- GLM 反事实预测超出 [0,1]：按原规则不截断，只要标准化均值有限且为正就保留点，同时输出越界诊断。
- 源文件在读取前后校验值不同：丢弃本次派生结果并标记 `source_changed`。

## 输出

```text
cells/<原 score_unit cell>/by_seed/stc_seed_000001.rds
report/combined_results.csv
report/paired_results.csv
report/method_summary.csv
report/paired_common_summary.csv
report/paired_differences.csv
report/input_audit.csv
report/missing_or_failed.csv
report/status_counts.csv
report/run_metadata.rds
configuration.rds
```

- `combined_results.csv`：STC 的 log RR、RR、SE、区间、风险平均、方差分量、状态、来源文件和 MD5。
- `paired_results.csv`：逐 seed 对应原 `score_unit` MAIC Wald 结果，保留原 MAIC 的不确定性说明。原 MAIC 区间不重新计算，两者并不自动具有相同的不确定性范围。
- `method_summary.csv`：方法各自成功样本上的偏倚、经验 SD、RMSE、覆盖率、宽度及返回率；返回且覆盖比例使用原计划分母。
- `paired_common_summary.csv` / `paired_differences.csv`：共同成功的配对比较，同时报告共同成功数量及原计划分母。`mean_squared_error_difference` 为 STC 减 MAIC 的均方误差差值。
- MAIC 基线会在两种 STC 区间分组中重复，以便对应比较；不是两套新 MAIC 分析。
- `input_audit.csv`：每个请求来源位置的读取状态、原始/处理后 MD5、完整标记和返回数量。
- `missing_or_failed.csv`：输入缺失、模型失败、目标支持不足和区间失败的明细，不会自动补 simulation。

`risk_C/risk_A`、`estimate/low/up` 等字段不改变原结果文件。报告中的目标真值仅用于计算误差/覆盖率，不参与拟合或 STC 标准化。输入全部缺失或没有保存真值时，误差和覆盖率为 NA，不从 DGP 补算真值。

## 单个 RDS 的调用

```r
source("stc_reuse.R")
original <- readRDS("/path/to/score_unit/by_seed/seed_000017.rds")
answer <- stc_from_record(original)
answer$summary       # 24 行 STC 结果
answer$paired        # 与原 MAIC 配对
answer$diagnostics   # 目标质量和、参数个数等
```

本地维护测试见源码目录 `tests/test_stc_reuse.R`，它使用本工作区原来保存的 validation RDS，不生成新数据；上传运行包不携带这些 RDS 或本地测试输出。

上传包仅包含两个运行源文件、本文档、验证记录、原工程的 BRM 上游许可说明及 SHA256 清单；不需要本地测试目录。保留 `BRMPLUS_LICENSE.txt` 以说明原工程相关公式/代码适配的上游来源。
