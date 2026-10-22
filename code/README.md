# RR / RD code corresponding to result_v2

当前绘图结果见 [../results](../results)。本目录仅保存源码、配置和代码说明。

本目录保存可共享的 RR、RD 完整源码，对应 `anchored_v1/result_cluster/result_v2` 的模型与算法。
整理日期：2026-10-22。历史结果、旧图、旧快照及补算归档未包含在内。

| 目录 | 内容 |
|---|---|
| [RR/code](RR/code) | 完整 RR v2 生成、估计、CI、模拟和汇总代码 |
| [RD/code](RD/code) | 完整 RD 代码及所需共享依赖 |
| [CI_methods](CI_methods/README_zh.md) | 单列非 Wald 区间：RR BC、RD BC、RD MN-ESS；含函数索引和源码审阅副本 |
| [plotting](plotting/README_zh.txt) | 绘图源码、单情景 CSV 与批量绘图入口 |
| [config](config/SCENARIOS.csv) | RR、RD 各 288 个情景的参数，不含结果 |
| [metadata](metadata) | 源码来源、版本与验证范围 |

详细说明见 [README_zh.txt](README_zh.txt)。RR 是 v2 DGM 加此前已发布 BC/FC 修复的整理副本；RD 保留 128 个冻结源码文件。
Wald 实现在完整主程序中，独立 CI 清单只列非 Wald 方法。RR 与 RD 的 BC 均采用参数 bootstrap 与 Blaker 型反演；RD 的额外文件还包括 MN-ESS 及模型/求解器适配代码。

## Quick checks

Run inside this directory:

```sh
python verify_checksums.py
Rscript --vanilla verify_runtime.R RR
Rscript --vanilla verify_runtime.R RD
Rscript --vanilla example_one_seed.R RR
Rscript --vanilla example_one_seed.R RD
```

The examples run point estimation for one seed with BC disabled. R and MASS are required; the RR FC rescue path also needs nleqslv. Full simulation and bootstrap CI runs were not repeated when preparing this bundle.

```sh
Rscript --vanilla plotting/plot_one_csv.R RR /path/to/combined_results.csv /path/to/new_RR_plots 1000
Rscript --vanilla plotting/plot_one_csv.R RD /path/to/anchored_rd_estimates.csv /path/to/new_RD_plots 1000
```

Plotting requires saved per-replicate CSV/RDS inputs, which are not included. Checksum verification covers the source sharing bundle; this GitHub README and the local `.gitattributes` are additional navigation/transport files.

Existing upstream license notices are preserved. They do not assign one new open-source license to all project code.
