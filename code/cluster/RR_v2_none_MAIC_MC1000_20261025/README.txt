RR v2 none+MAIC -- eight scenarios, 1000 data seeds each -- 2026-10-25

运行设定
--------
不在结局回归加入 score；用 MAIC 权重拟合结局模型：
MAIC_EST_MODE=none，MAIC_MODEL_WEIGHTING=maic。
MAIC 校准 all，即 z1..z6；diffuse；high EM；当前 v2 数据生成机制。
原正式程序的全部方法、方差分支和 P 值计算保留。
CI 粗搜索/确认/精修=600/1000/1000，P 值 B=500。
保留最终版多起点优化、nlminb 备用求解和四处浮点阈值修正。
边界截断是有效的方法结果，不作为失败重算。

task  event    n    名义 Full ESS
1     common  100  .15
2     common  100  .50
3     common  500  .15
4     common  500  .50
5     rare    100  .15
6     rare    100  .50
7     rare    500  .15
8     rare    500  .50
n 指每个试验的总人数，AC/BC 各 n 人，每个观察臂 n/2 人。
每种情形 data seed=1..1000；CI seed=1000000+data seed。
同一 v2 情形、相同 data seed 的原始数据生成过程与此前相同。
本包重新生成并进行 MAIC 加权分析，不复用旧 none+unit 的估计结果。

资源与目录
----------
每组一个 Slurm 数组元素，完整处理1000 seed；64核、10G、4h。
不分块、不设置数组并发上限；实际启动仍由调度器安排。
代码目录：$HOME/RR_v2_none_MAIC_MC1000_20261025
独立结果目录：${SCRATCH:-/scratch/$USER}/RR_v2_none_MAIC_MC1000_results_20261025
任务名称 RRnoneMAIC_g1；只含这8组。
使用独立 RRNMW_* 环境变量，不继承旧 RRFULL_/RRSTAB_/RDSTAB_ 的结果目录。
默认账户 def-liteep，可 export RRNMW_ACCOUNT 更改。
本包包含完整代码，不依赖 Nibi 上已安装的旧代码包或 Rorqual 输入。
113个运行源文件与原 v2 全方法包完全一致；只改变分析清单和独立运行入口。

本机 PowerShell 上传
-------------------
$pkg = "D:\codexresults\ITC\rorqual\compareITC\anchored\RR_v2_none_MAIC_MC1000_20261025.tar.gz"
scp -o ConnectTimeout=60 -o ConnectionAttempts=3 "$pkg" "${pkg}.sha256" yanjin41@nibi.alliancecan.ca:~/
if ($LASTEXITCODE -ne 0) { throw "none+MAIC代码上传失败" }
ssh yanjin41@nibi.alliancecan.ca

Nibi 提交（含前面要求的自动检查续跑）
-----------------------------------
cd "$HOME" &&
sha256sum -c RR_v2_none_MAIC_MC1000_20261025.tar.gz.sha256 &&
tar -xzf RR_v2_none_MAIC_MC1000_20261025.tar.gz &&
cd "$HOME/RR_v2_none_MAIC_MC1000_20261025" &&
bash submit_auto.sh

这会先提交一个8元素计算数组、一个 afterany 汇总任务（1核8G1h），
再挂一个 afterany 检查器（1核4G30min）。无需保持SSH连接。
若只需要一轮计算，改用 bash submit.sh；不要同时运行两条初次提交命令。

自动续跑也只处理本包的 none+MAIC，64核/10G/4h保持不变：
完整 seed 直接复用；缺失 seed 补算；中断 seed 先保留 partial_history，
再从该 seed 的起点重算。不会从某个 bootstrap 循环的中间位置接着算。
完整处理但某方法区间失败的 seed 不自动重抽、不反复重算。
默认最多补10轮；连续两轮没有新增完整 seed 或发现损坏/不兼容RDS时停止。
自动续跑结束后会排队刷新完整汇总。检查/补跑ID记在本独立结果目录中。
进度入口：auto_none_MAIC_20261025/latest_chain.txt。

手动续跑、只补未完成部分
-----------------------
bash "$HOME/RR_v2_none_MAIC_MC1000_20261025/submit.sh"
已有本包作业运行时拒绝重复提交。原 none+unit 与其旧自动链不受此包控制。
若已手动提交本包，后来要挂自动续跑，可执行：
bash "$HOME/RR_v2_none_MAIC_MC1000_20261025/auto/submit_auto_none.sh"
只补指定情形时用 submit.sh 1,2；自动检查入口始终覆盖全部8组。

结果与验证
----------
各组 cells/<cell>/by_seed/seed_000001.rds 等边跑边保存；complete=TRUE才是完整处理。
report/progress.csv：逐情形完整/中断/缺失记录。
report/combined_results.csv：所有方法的完整结果。
report/all_methods_summary.csv：方法返回、有效区间覆盖率、固定1000分母的返回且覆盖比例等。
边际结果尺度为 log RR，目标为 BC 总体 A/C。
VALIDATION.txt 记录实际完成的本地检查。数值模型仅用CI替身检查加权接线，
不另跑coverage先导试验，不声称已完成生产的8000个区间。
当前仅本地准备，没有上传、提交、修改或取消任何 Nibi 作业。
