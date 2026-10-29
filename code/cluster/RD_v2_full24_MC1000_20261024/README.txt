RD v2 正式运行：8 个数据情形 x 3 种分析设定 x 1000 seed
版本：RD_v2_full24_MC1000_20261024
状态：已准备本地包；未上传或提交 Nibi。

范围与来源
当前 RD_final_B600_F1000_20261024/code 的 147 个运行源码逐字节复制，不再修改统计代码。
不使用历史 RD_effect_v2 的旧区间 helper，也不读取冻结 bootstrap 补测的输入/输出。
从当前 RD v2 DGM 重新生成模拟数据，运行 RD 原正式程序全部九种方法及 P 值。
无需从 Rorqual 搬输入、seed 或旧结果，只需上传此包及校验文件。

情形顺序（n 是每个试验的总样本量）
1 Common n=100 ESS=.15；2 Common n=100 ESS=.50；
3 Common n=500 ESS=.15；4 Common n=500 ESS=.50；
5 Rare n=100 ESS=.15；6 Rare n=100 ESS=.50；
7 Rare n=500 ESS=.15；8 Rare n=500 ESS=.50。
AC 和 BC 各有 n 人；每个试验两臂均分，n=100 时每臂 50 人。
全部 diffuse、high EM、MAIC 校准 all（z1..z6）。
ESS 是设计的名义 Full ESS 比例；实际样本 ESS 保存在每个 seed 的 estimates 中。
当前 RD v2 以加性风险直接生成：BC 人群的 A-B RD 为 common=.14、rare=.05。
BC C/A/B 风险为：common .418596226792393/.642197586315061/.502197586315061；
rare .0488309131792698/.1089376565366341/.0589376565366341。
强制检查 dgm=webster2024_rd_ess_em_v2、scenario version=2。
RD 和 RR 使用各自的 DGM 与效应尺度；不能将 RR 的风险参数直接代入 RD。

每个情形三种分析（tasks.csv 中按此顺序）
none：estimation_mode=none，model_weighting=unit，结局模型不加 score。
score_unit：estimation_mode=scores，model_weighting=unit。
score_MAIC：estimation_mode=scores，model_weighting=maic。
三组的 MAIC 标准化和相应 MAIC 方法均使用 all 校准；none 仍计算 MAIC 权重。
score 用 z1,z2,z3 作为 EM，z4,z5,z6 作为 non-EM，bc_estimated 标准化。
brm 使用当前正式的 linear h_E on atanh(RD) 模型；none 下为 intercept only。
同一情形三组用相同 data seed=1..1000，生成数据一致；CI seed=1000000+data seed。
共 24000 个分析记录，等于 8000 份生成数据分别在三组中分析。
九种原方法：brm、brm-BC、brm-FC、brm-FC-BC、IB、RIB、RLPM、MN_eff、RD_MAIC_sandwich。

数值设定
CI 粗搜索/确认/精修=600/1000/1000；独立 P 值预算=500。
四处判定统一 acceptability >= (1-level)-1e-12。
保留 L-BFGS-B、两次续算、最多 3 个确定性不同初值、nlminb 和 RD 稳定边界恢复。
保留相同数据、目标函数、约束及接受规则，不丢弃或重抽失败的 bootstrap 响应。
每次 exact 区间调用的内部 time.limit.sec=10800（3h），独立于整个作业的 Slurm 时长。
这项数值上限沿用最终版；更长 Slurm 时长不自动解除单次求区间上限。

资源配置与证据
none：n=100/500 均 32 核、4G、8 小时，按用户 2026-10-25 要求调整。
score_unit/score_MAIC：n=100 为 96 核、20G、18 小时；n=500 为 96 核、20G、19 小时。
每个情形/分析一个完整 1000 seed 数组元素，不分 seed 块，不设置 % 并发限制。
默认按三种资源组合形成 3 个数组，每个数组 8 个元素，共 24 个计算任务；
另有 1 核、8G、1 小时 afterany 汇总任务。实际何时启动由 Nibi 调度器决定。
--mem=4G/20G 表示单节点整个任务的内存，不是每核内存。

按用户要求，score 配置依据 RR 的既有记录选择，不声称是 RD 完整流程实测。
resource_evidence_B.csv 原样保留 6 个 RR B 档作业记录；其中 200 seed/64 核耗时
1h40m--2h31m，.batch MaxRSS 为 4.3--6.9 GiB。
折算到 1000 seed/96 核，单 FC 分支约 5.6--8.4h。
完整 score 程序有 brm-BC 和 brm-FC-BC 两个主要 Exact 区间分支，还需生成数据、
拟合及计算 P 值。RR 配置按两个分支折算，再留约 2h，得到 n100=18h、n500=19h。
这些旧 RR 情形并非全部匹配新 diffuse/high/all/ESS 条件，score_MAIC 也无直接计时。
RD 的数值边界行为和权重可能另增耗时，因此上述数值是申请上限与规划估计，
不能保证 1000 seed 一次完成；程序完成立即退出，超时后可按同版本续跑。
如果需要进一步收紧/提高 RD 时长，请提供旧 RD B 档作业的 Elapsed、AllocCPUS、
State、Timelimit 和 .batch MaxRSS；不需要为此另做 bootstrap pilot。
Slurm 参数参考：https://slurm.schedmd.com/sbatch.html

本机 PowerShell：上传并登录
$pkg = "D:\codexresults\ITC\rorqual\compareITC\anchored\RD_v2_full24_MC1000_20261024.tar.gz"
scp -o ConnectTimeout=60 -o ConnectionAttempts=3 "$pkg" "${pkg}.sha256" yanjin41@nibi.alliancecan.ca:~/
if ($LASTEXITCODE -ne 0) { throw "RD正式版上传失败" }
ssh yanjin41@nibi.alliancecan.ca

Nibi：校验、解压和一次提交全部 24 组
cd "$HOME" &&
sha256sum -c RD_v2_full24_MC1000_20261024.tar.gz.sha256 &&
tar -xzf RD_v2_full24_MC1000_20261024.tar.gz &&
cd "$HOME/RD_v2_full24_MC1000_20261024" &&
bash submit.sh

账户默认 def-liteep，可用 RDFULL_ACCOUNT 更改。
2026-10-25 已排队的 none 数组为 24958765（原 16 核、4h）。在 Nibi 直接修改：
scontrol update JobId=24958765 NumCPUs=32 CPUsPerTask=32 MinCPUsNode=32 TimeLimit=08:00:00
随后核对：squeue -j 24958765 -o "%.24i %.18j %.10T %.6C %.10l %.8m %R"
以上命令由用户在 Nibi 执行；本地打包过程未连接集群或修改实际任务。
无需为已排队作业取消或重交。若服务器仍是旧包，之后续跑时显式使用：
RDFULL_NONE_CPUS=32 RDFULL_NONE_TIME=08:00:00 bash "$HOME/RD_v2_full24_MC1000_20261024/submit.sh"
输出默认 ${SCRATCH:-/scratch/$USER}/RD_v2_full24_MC1000_results_20261024。
每个任务用一个节点，--ntasks=1，R fork 并行独立 seed；BLAS/OMP 每个 worker 单线程。
加载模块 StdEnv/2023 gcc/12.3 r/4.3.1。
可设置 RDFULL_NONE_{CPUS,MEM,TIME}、RDFULL_SCORE_{CPUS,MEM,TIME_N100,TIME_N500}、
RDFULL_MAIC_{CPUS,MEM,TIME_N100,TIME_N500} 覆盖资源；不影响模拟及断点版本。

断点续跑
再次 bash "$HOME/RD_v2_full24_MC1000_20261024/submit.sh"。
已有完整组会跳过；剩余组复用同版本处理完成的 seed。
活动作业存在时拒绝重复提交，避免争用同一结果目录。
部分 seed 的旧检查点先逐字节保存在 partial_history，然后从该 seed 的原始起点重算。
完整模拟中某个方法/区间失败的记录保留，不自动删除、重抽或反复求到成功。
要运行全新一轮请指定新的 RDFULL_OUT。
可提交指定 task_id，例如 bash submit.sh 1,2,3 仅对应情形 1 的三种分析。
task_id=3*(scenario_id-1)+分析序号（none=1、score_unit=2、score_MAIC=3）。

结果
cells/<cell>/by_seed/rd_seed_00000001.rds 等为完整逐 seed 检查点。
各组保留 anchored_rd.rds、anchored_rd_estimates.csv、anchored_rd_coefficients.csv、
anchored_rd_summary.csv、原方法状态表及 missing_seeds.csv、incomplete_seeds.csv。
report/progress.csv 记录每组已处理、部分记录及缺失 seed 数；未完成汇总任务返回非零。
report/combined_results.csv 合并九种方法并附 task_id/scenario_id/analysis。
report/all_methods_summary.csv 提供各方法有效区间数、失败数、coverage、MCSE、宽度、bias。
coverage_valid_pct 分母为有效区间数；returned_and_covered_pct 分母固定 1000。
BC_boundary_truncated 作为有效区间保留。CASE_PROCESSED 只表示 1000 条模拟已处理，
并不表示每个方法都成功返回区间。VERIFY.txt 区分接口验证与实际正式模拟。
