# 当前 Nibi 结果对应代码

RR 的 113 个运行文件、RD 的 147 个运行文件，均与六个当前 CSV 导出包的 `RUNTIME_FILES.csv` 中 SHA256 完全一致。
STC 的两个程序 MD5 与保存的运行配置一致。证据见 [metadata/RUNTIME_VERIFICATION.json](metadata/RUNTIME_VERIFICATION.json) 和 [metadata/SOURCE_FILES.csv](metadata/SOURCE_FILES.csv)。

| 功能 | RR | RD |
|---|---|---|
| 估计流程 | [cluster_all_methods.R](RR/code/estimation/cluster_all_methods.R) → [anchored_comparison.R](RR/code/estimation/anchored_comparison.R) | [rd_cluster.R](RD/code/RD/estimation/rd_cluster.R) → [rd_comparison.R](RD/code/RD/estimation/rd_comparison.R) |
| BC 区间核心 | [CI_marginal_target.R](RR/code/estimation/helpers/CI_marginal_target.R) | [CI_marginal_rd_target.R](RD/code/RD/estimation/CI_marginal_rd_target.R) |
| Coverage 汇总 | [plot_metrics()](plotting/RR/plot_cluster_results.R) | [rdplot_metrics()](plotting/RD/plot_rd_cluster_results.R) |

- RD 的 MN-ESS： [rd_methods.R](RD/code/RD/estimation/rd_methods.R)。
- 非 Wald CI 的独立审阅入口：[CI_methods](CI_methods/README_zh.md)。
- STC 运行入口：[run_stc_reuse.R](STC_RR/run_stc_reuse.R)，核心实现：[stc_reuse.R](STC_RR/stc_reuse.R)。
- 统一重绘入口：[plot_nibi_results.R](plotting/plot_nibi_results.R)。
- [config/SCENARIOS.csv](config/SCENARIOS.csv)：当前 48 个 RR/RD 分析情景。
- [cluster](cluster)：完整、保持原字节的 full24 和 none_MAIC 运行包；各包含自身代码、配置、运行入口与原 SHA256SUMS。运行前阅读包内 README；不会在本地自动提交任何作业。

当前 BC 默认参数为 coarse / confirm / refine / p-value = 600 / 1000 / 1000 / 500，置信水平 .95，二分阈值 .001，二分上限 15。
完整 full24 原包包含 `none`（单位拟合权重）；当前主结果采用单独的 `none_MAIC` 导出，两者不得混用。

源码检查（在本目录运行）：

```sh
python verify_checksums.py
Rscript --vanilla verify_runtime.R RR
Rscript --vanilla verify_runtime.R RD
```

从仓库根目录重绘：

```sh
Rscript --vanilla code/plotting/plot_nibi_results.R .
```

运行检查只解析和加载源码、核对 DGM 与区间默认参数，不运行 MC1000 或 bootstrap。
`RR/code`、`RD/code` 与 cluster 包内的数值源码保持运行快照原字节；整理说明、校验脚本和绘图驱动为新增共享入口。
原有 BRMPLUS 等许可声明均保留。
