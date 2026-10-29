最新 Nibi 结果重绘

在仓库根目录执行：
Rscript --vanilla code/plotting/plot_nibi_results.R .

输入：original result/config/SCENARIOS.csv 指定的48份逐次CSV，及STC_RR/report/combined_results.csv.gz。
RR、RD各三种分析，每种8个情景。STC的fixed_bc_summary和bc_mean_delta分别绘制8个情景。
输出：每个情景bias.png、se.png、accuracy.png、coverage.png、四页metrics.pdf、metrics.csv和analysis_notes.txt。
RR/RD汇总表为results/metrics_all.csv，STC为results/STC_RR/metrics_all_runs.csv。

不生成数据、不拟合模型、不计算新CI。使用保存的区间端点与真值重新汇总。
RR在log(RR)尺度，RD在自然RD尺度。Coverage以有效返回且有真值的区间为分母。
Accuracy=平均有效模型SE/同一有效SE子集上的估计经验SD。
原始失败行保留；没有模型SE的BC/MN-ESS不进入Accuracy。
RR图延续既有展示规则，省略brmad、LP-BFC、RLP-BFC；指标表和原始CSV保留所有方法。
RR Bias/SE显示范围默认5/4，RD默认1/1；图面裁剪不改变计算结果。

每次绘图会将RR/RD核心指标与独立CSV核对表比较，STC与其源method_summary.csv比较。
validation目录中的独立参考表与当前输入包SHA一致；不重新运行模型验证。
