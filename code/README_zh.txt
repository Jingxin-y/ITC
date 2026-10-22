RR / RD 源码分享包 — result_v2 对应版本
2026-10-22 修订

目录
RR/code：完整 RR v2 数据生成、MAIC、估计、置信区间、模拟和汇总代码。
RD/code：完整 RD 源码及其共享依赖；RD 专用模块在 RD/code/RD。
CI_methods：单独列出的非 Wald 区间说明、函数索引和原样源码审阅副本。
  包括 RR BC、RD BC，以及 RD MN-ESS；不单列 Wald。
plotting：独立绘图源码、单情景 CSV 绘图入口、批量绘图入口及输入说明。
config/SCENARIOS.csv：RR、RD 各 288 个情景的参数配置，不含估计结果。
metadata：本分享包的源码来源、文件 SHA256 和检查摘要。

此修订不包含历史结果表、旧图、审计报告、旧源码快照、修复案例输入或补算归档。

版本
RR 保留 v2 生成机制，使用此前已发布的 BC 与 FC 求解器修复。
它是按对应源文件整理的可运行副本，并非经认证的单一集群最终快照。
RD 的 128 个冻结源码文件与既有 SHA256 清单一致，另附许可文件。
全部主程序源码与上一分享包相同；本次修改分享结构、说明和绘图入口。
此包对应 result_cluster/result_v2 的模型/算法，不包含 full_CI_300_500 试验算法。

先阅读
1. CI_methods/README_zh.md：非 Wald 区间的目标、计算步骤、参数和检查入口。
2. CI_methods/METHOD_INDEX.csv：直接定位函数、源码文件和行号。
3. plotting/README_zh.txt：如何从保存的 CSV/RDS 重绘四类图。

快速检查（在分享包根目录运行）
python verify_checksums.py
Rscript --vanilla verify_runtime.R RR
Rscript --vanilla verify_runtime.R RD
Rscript --vanilla example_one_seed.R RR
Rscript --vanilla example_one_seed.R RD
后两条仅为一个 seed 的点估计接口示例，不计算 BC，也不运行 MC1000。
主要依赖 R 和 MASS；RR FC 后备求根需 nleqslv（原验证版本 3.3.5）。

绘图示例
Rscript --vanilla plotting/plot_one_csv.R RR /path/to/combined_results.csv /path/to/new_RR_plots 1000
Rscript --vanilla plotting/plot_one_csv.R RD /path/to/anchored_rd_estimates.csv /path/to/new_RD_plots 1000
批量入口和完整输入要求见 plotting/README_zh.txt。

正式配置
每尺度 288 个配置：scores/unit 144 个、none/MAIC 144 个；每配置 seeds=1:1000。
RR v2 必须使用包内 RR 生成器；它与本地 anchored_v1 根目录的 RR v1 不同。
RD scores/unit 使用 linear_atanh；none/MAIC 的有效模型为 intercept_only。
源码 simulation_config 中的 seeds=1:10 是演示设置，正式参数见 config/SCENARIOS.csv。

范围
包内不附模拟结果。绘图需提供原始逐次估计 CSV 或 RDS，只有最终 metrics.csv
不能恢复 Bias/SE 箱线图的逐次分布。本次没有重新运行全部模拟或完整 BC。
保留 BRMPLUS_LICENSE.txt 等已有声明；该声明不代表全部项目代码统一获开源许可。
