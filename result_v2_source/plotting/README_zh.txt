绘图代码与运行入口

plot_one_csv.R：可迁移的单情景入口，读取保存的逐次估计 CSV。
RR/plot_cluster_results.R：这批 RR 图形使用的原始 helper，含指标与绘图函数。
RD/plot_rd_cluster_results.R：这批 RD 图形使用的冻结 helper。
plot_runs.R：调用原 helper，对完整结果目录批量绘图。
所有入口只读取结果，不重新生成数据、拟合模型或计算置信区间。

单情景 CSV（在包根目录运行）
Rscript --vanilla plotting/plot_one_csv.R RR /path/to/combined_results.csv /path/to/new_RR_plots 1000
Rscript --vanilla plotting/plot_one_csv.R RD /path/to/anchored_rd_estimates.csv /path/to/new_RD_plots 1000
最后参数为计划重复次数，默认1000；不是从成功记录数推定分母。
输入必须是一个情景的逐 seed 长表，每个 seed 保留所有方法/区间记录，
失败记录也应保留其状态。输入当前采用 RR 每seed18条、RD每seed9条的原始导出格式。
列要求和情景一致性由入口检查；多个情景不能拼在同一 CSV 中。

批量结果目录
Rscript --vanilla plotting/plot_runs.R RR /path/to/RR_results --out=/path/to/new_RR_plots
Rscript --vanilla plotting/plot_runs.R RD /path/to/RD_results --out=/path/to/new_RD_plots
可追加 --batch=文本、--format=png,pdf 或 --dry-run。
RR 使用有 by_seed/ 的运行目录，可读取逐 seed RDS 或既有汇总CSV；
RD 使用包含 anchored_rd.rds 的单情景目录，或其父目录。
迁移时直接传入实际输入输出目录，无需原集群上的 HOME/scratch 指针。

每情景输出
bias.png、se.png、accuracy.png、coverage.png；metrics.pdf（四页）；metrics.csv。
Bias 与 SE 使用逐次估计的箱线图；Accuracy 与 Coverage 使用指标图。

统计和显示口径沿用原 helper
RR 在 log(RR) 尺度；RD 在自然 RD 尺度。
Accuracy = 平均 SE / 估计量的经验SD，分子分母使用相同的SE有效重复。
Coverage 以有效返回且有对应真值的区间为分母，同时报告区间返回比例。
BC 和 MN-ESS 没有模型 SE，因此不进入 Accuracy 图。
RR 继续排除 brmad、LP-BFC、RLP-BFC 的展示；其计算代码仍在完整主程序中。
重复的点估计在箱线图中仅展示一次；边界截断区间保留原状态。
显示上限仅裁剪图面，不改变均值、SD、RMSE或覆盖率的计算。
RR Bias/SE 显示临界值默认5/4；RD 默认1/1。

输入说明
最终汇总 metrics.csv 只有统计量，不能重建 Bias/SE 箱线图。
需要保存的逐次估计 CSV 或 RDS；本分享包不附历史结果或测试生成的图片。
新入口复用原 helper 的指标与绘制函数；新增部分仅做参数、输入检查和字段衔接。
