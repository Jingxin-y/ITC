RR v2 正式运行：8 个数据情形 x 3 种分析设定 x 1000 seed
版本 RR_v2_full24_MC1000_20261024

Material Passport
任务：按用户指定设计准备 RR 正式 Monte Carlo 程序。
状态：PREPARED；尚未上传或提交 Nibi；资源根据既有 RR 作业线性外推，并按用户最新偏好收紧。
来源：本地 RR_effect_v2_20261015 当前源码，包含已同步的 RR 最终区间算法。
原始数据/历史结果：不改动；此包从 v2 机制重新生成数据。

设计
----
scenarios.csv 严格按用户次序列出 8 种情形：
1 Common n=100 ESS=.15；2 Common n=100 ESS=.50；
3 Common n=500 ESS=.15；4 Common n=500 ESS=.50；
5 Rare   n=100 ESS=.15；6 Rare   n=100 ESS=.50；
7 Rare   n=500 ESS=.15；8 Rare   n=500 ESS=.50。
每个试验的总人数为 n，两个观察臂各 n/2；AC 和 BC 各有 n 人。
全部使用 v2、diffuse、high EM、MAIC 校准 all（z1..z6）。
ESS 是名义 Full ESS 比例，模拟样本的实际 ESS 保留在输出中。

每种情形的三组分析依次为：
none       ：MAIC_EST_MODE=none；MAIC_MODEL_WEIGHTING=unit。
score_unit ：MAIC_EST_MODE=scores；MAIC_MODEL_WEIGHTING=unit。
score_MAIC ：MAIC_EST_MODE=scores；MAIC_MODEL_WEIGHTING=maic。
none 不在结局模型加入 score；MAIC 的 all 校准仍用于标准化及 CMH。
score 的 EM 为 z1,z2,z3，non-EM 为 z4,z5,z6。
同情形三组使用相同 data seed=1..1000，数据完全相同；CI seed=1000000+data seed。
共 24000 个分析记录，对应 8000 份生成数据在三种分析下的比较。
运行原正式程序的全部方法，保留正式模式的 P 值计算。

数值与 DGM 版本
---------------
CI 粗搜索/确认/精修=600/1000/1000，P 值 nsim=500。
先原 L-BFGS-B 和两次续算，再最多 3 个确定性不同初值，最后 nlminb。
四处阈值为 acceptability >= (1-level)-1e-12。
CI helper 与 RR_final_B600_F1000_20261024 的最终 helper 字节完全相同。
所有 data/estimation/summary 源码均来自当前 RR_effect_v2_20261015，未再修改。
v2 的 BC 风险 C/A/B：
common .418596226792393/.514669087904676/.60；
rare .0488309131792698/.10/.15；baseline_rr_E=.59。
强制检查 dgm=webster2024_rr_ess_em_v2、scenario version=2。
不读取旧 frozen input，也不复制或复用旧模拟结果；无需从 Rorqual 搬数据。

资源建议
--------
每组一个 Slurm 数组元素，完整处理 1000 seed，不分块，不设数组并发上限。
none：所有样本量均 16 核、4G、4 小时，按用户指定申请。
score_unit：n=100 为 96 核、20G、18 小时；n=500 为 96 核、20G、19 小时。
score_MAIC：先用与 score_unit 相同资源，无直接计时证据，不另加大余量。
默认分成三个资源数组，各 8 个元素，共 24 个计算任务。
另有一个 afterany 汇总任务：1 核、8G、1 小时。
这些是申请上限，不是预计必然用满的时间；程序完成即退出。
根据用户提供的已运行 RR B 档，而非 RD：200 seed/64 核约 1h40m--2h31m，
记录的峰值内存约 4.3--6.9 GiB；按线性吞吐外推，1000 seed/96 核的
单个 FC 区间分支约 5.6--8.4h，仅用于粗略参照。
正式程序还有 brm、brm-FC，none 下可能还有 brmad 的 Exact 分支及 P 值，
并重新生成数据和拟合所有模型。新 ESS/high/all、MAIC 拟合和多起点恢复
都会影响实际耗时，因此不能把旧补测时间当作正式时长的精确预测。
score_unit 的折算明确计入两个主要 Exact 分支（brm 与 brm-FC）：
n=100 用原 n=100/200 的 B 档较慢记录约 2h20m，乘 1000/200、乘 64/96、
再乘两个主要分支，约 15.6h；增加约 2h 后向上取整为 18h。
n=500 用 B 档 2h31m，同样折算约 16.8h；增加约 2h 后向上取整为 19h。
两个分支实际耗时未必相同，上式是透明的规划近似，不是完整流程实测。
none 按用户经验指定 16 核/4h；完整 none/MAIC 用时和内存尚未实测，超时可续跑。
4/20G 是资源配置值，未声称它们是实测最低内存。
resource_evidence_B.csv 保存采用的 6 个 RR B 作业记录及其原情形。
MaxRSS 来自作业 .batch 行，空白作业主行不能当成零内存。
参考 Slurm 原始说明：https://slurm.schedmd.com/sacct.html

本机 PowerShell
---------------
$pkg = "D:\codexresults\ITC\rorqual\compareITC\anchored\RR_v2_full24_MC1000_20261024.tar.gz"
scp -o ConnectTimeout=60 -o ConnectionAttempts=3 "$pkg" "${pkg}.sha256" yanjin41@nibi.alliancecan.ca:~/
if ($LASTEXITCODE -ne 0) { throw "RR正式版上传失败" }
ssh yanjin41@nibi.alliancecan.ca

Nibi 提交
---------
cd "$HOME" &&
sha256sum -c RR_v2_full24_MC1000_20261024.tar.gz.sha256 &&
tar -xzf RR_v2_full24_MC1000_20261024.tar.gz &&
cd "$HOME/RR_v2_full24_MC1000_20261024" &&
bash submit.sh

默认输出：${SCRATCH:-/scratch/$USER}/RR_v2_full24_MC1000_results_20261024
账户默认 def-liteep，可用 RRFULL_ACCOUNT 更改。
资源默认值见 resources.csv；可分别覆盖：
none：RRFULL_NONE_CPUS、RRFULL_NONE_MEM、RRFULL_NONE_TIME。
score_unit：RRFULL_SCORE_CPUS、RRFULL_SCORE_MEM、RRFULL_SCORE_TIME_N100/N500。
score_MAIC：RRFULL_MAIC_CPUS、RRFULL_MAIC_MEM、RRFULL_MAIC_TIME_N100/N500。
例如 RRFULL_NONE_TIME=04:00:00；walltime 使用 Slurm 标准格式。
覆盖值只影响资源，不改变模拟设定。不要沿用旧 MAIC/RDSTAB/RRSTAB 的提交入口。

后续续跑
--------
bash "$HOME/RR_v2_full24_MC1000_20261024/submit.sh"
会跳过已有处理完成标记的整组；其余组重新提交，组内复用同版本的 complete seed。
有活动作业时拒绝重复提交，以免争用同一输出目录。
未完成的部分 seed 先保留 partial_history，再按原算法从该 seed 起点重算。
已完成但方法/区间失败的记录保留，不自动重抽或反复计算。
如需要全新的整轮重跑，请显式设置一个新的 RRFULL_OUT 目录。
可只续跑指定任务，例如情形 1 的三组：bash submit.sh 1,2,3。
task_id=3*(scenario_id-1)+分析序号（none=1，score_unit=2，score_MAIC=3）。

结果
----
各组 cells/<cell>/by_seed/seed_000001.rds 等为完整逐 seed 检查点。
原正式方法的 combined_results.csv、combined_coefficients.csv、缺失/未完成清单保留。
report/progress.csv 区分已处理、未完成和缺失 seed。
report/combined_results.csv 合并所有组，并附 task_id/scenario_id/analysis。
report/all_methods_summary.csv 按方法给出返回数、coverage、Monte Carlo SE、宽度与 bias。
coverage_valid_pct 的分母为有效区间数；returned_and_covered_pct 的分母固定 1000。
CASE_PROCESSED 表示 1000 条模拟记录处理完毕，不代表所有方法均返回有效区间。
边界截断是方法正常结果，作为有效区间保留。诊断失败不删除 bootstrap 响应。
VERIFY.txt 记录本地实际完成的配置/代码检查；不包含新的 coverage 先导试验。
