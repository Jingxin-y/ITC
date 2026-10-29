RR brm-FC 失败后备求根（2026-10-19）

现有 source_anchored() 会自动加载本实现，无需另外 source 独立补丁。
实现位置：estimation/helpers/MLE_Point_Firth_for_RR.R。
只有 RR 入口新增失败后备；共享 FC 主算法、RR 得分/导数、概率映射和 RD 保持原样。

执行顺序
1. 先完整执行原 max.likelihood.firth() 主算法。
2. 若主算法 convergence=TRUE，立即原样返回完整对象，完全不调用后备。
3. 若未收敛且 max.step-primary$step>0，从主算法保存的最佳有限候选继续：
   nleqslv(method="Newton", global="none", stepmax=5,
           ftol=min(thres,1e-10), xtol=1e-12, allowSingular=TRUE)。
   若 Newton 未通过且还剩预算，再从同一候选尝试局部 Broyden。
4. 最终必须独立通过原始二元响应、归一化拟合权重、原 score_scale 上的
   max(abs(adjusted.score/score_scale))<=thres 验收；两臂概率处于原
   (1e-8,1-1e-8) 内部；参数、得分、信息和协方差有限且矩阵正定。
   不依靠求解器返回码认定成功。仍失败则保留原失败候选和失败状态。

预算与诊断
max.step 仍然是主算法和后备合计的总迭代上限。原默认 300 和 thres=1e-6 保持不变。
当后备包不可用、没有剩余预算、或找不到合格根时，不把失败结果标记为收敛。
底层 result$diagnostics$rescue 记录尝试、方法、包版本、迭代/得分评估数和验收；
MLEst 返回对象中对应 fit$fc.diagnostics$rescue。
primary_diagnostics 保存未修改的主算法诊断；成功主流程的对象不新增任何字段。
后备没有随机初始化，不抽取随机数。

依赖
原本成功的拟合不需要新依赖。使用失败后备时，当前 R library 必须有 nleqslv：
  install.packages("nleqslv", repos="https://cloud.r-project.org")
可用 requireNamespace("nleqslv", quietly=TRUE) 检查当前 R 是否能找到它。
使用用户 library 即可；本实现不会主动安装包或修改 .libPaths()。
本地验证版本为 R 4.3.2 / nleqslv 3.3.5；集群应在实际 R 环境中复核。

复现
在 anchored_v1 目录运行：
  Rscript estimation/tests/check_rr_fc_rescue.R . path/to/failed_fc_details.rds output_dir
可选第四参数为合并前的 MLE_Point_Firth_for_RR.R；第五参数为历史 checkpoint
路径清单（UTF-8，每行一个 RDS）。不提供第四参数时，基准为未修改的共享主算法。
该测试同时检查：八条失败经正式 MLEst 入口收敛、原始方程/概率/协方差验收、
已收敛拟合完整返回对象和 RNG 一致、后备未被调用、短预算、缺包和权重缩放。

点估计及一阶协方差恢复，不代表已重跑 BC 区间或总体 bias/RMSE/coverage。
新增有效记录之后需要通过原分析流程重新汇总。原始模拟 checkpoint 不应被手工覆盖。
本次八条输入的成功不构成任意样本、任意初值下的全局收敛或根唯一性保证。
