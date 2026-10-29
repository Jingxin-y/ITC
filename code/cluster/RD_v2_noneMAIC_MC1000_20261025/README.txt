RD v2 none+MAIC 正式包：8 种情形，每种 1000 seed
版本：RD_v2_noneMAIC_MC1000_20261025
状态：本地准备；尚未上传或提交 Nibi。

分析定义
estimation_mode=none：结局回归不加入 score 或协变量，只保留处理效应和截距。
model_weighting=maic：brm、brm-FC、IB/RIB、RLPM 的拟合使用 MAIC 权重。
预测风险/效应仍用 MAIC 权重标准化到 BC 目标人群。
none 的恒定预测效应在标准化前后相同；此包的关键变化是拟合阶段已使用 MAIC 权重。
MN_eff 与 RD_MAIC_sandwich 同样使用 MAIC 权重计算两臂风险。
MAIC 校准使用全部六个协变量 z1..z6；BC 汇总均值用于计算权重。
brm 与 brm-FC 沿用当前 RD 的回归、点估计优化器、收敛阈值和 FC 修正，不替换统计方法。
保留原正式九种方法：brm、brm-BC、brm-FC、brm-FC-BC、IB、RIB、RLPM、MN_eff、RD_MAIC_sandwich。

情形与随机数
全部采用当前 RD v2（webster2024_rd_ess_em_v2）、diffuse、high EM、all。
task_id / scenario_id 顺序：
1 Common n=100 ESS=.15；2 Common n=100 ESS=.50；
3 Common n=500 ESS=.15；4 Common n=500 ESS=.50；
5 Rare n=100 ESS=.15；6 Rare n=100 ESS=.50；
7 Rare n=500 ESS=.15；8 Rare n=500 ESS=.50。
n_AC=n_BC=n 为每个试验总人数，每臂 n/2；n=100 时每个试验两臂各 50 人。
ESS 是名义 Full ESS 比例，实际样本 ESS 保存在输出中。
每种情形使用 data seed=1..1000，CI seed=1000000+data seed。
与原 full24 包相同情形/seed 的生成数据一致，便于配对比较。
此包独立生成数据，不读取或复用旧 none+unit 或 score 组的结果。

算法与来源
147 个 code/ 运行文件逐字节复制自 RD_v2_full24_MC1000_20261024/code，
也与 RD_final_B600_F1000_20261024/code 一致；科学源码没有修改。
正式 CI 粗搜索/确认/精修仍为 B=600/1000/1000；独立 P 值预算=500。
四处接受阈值统一为 acceptability >= (1-level)-1e-12。
保留 L-BFGS-B、续算、最多 3 个确定性替代初值、nlminb 与 RD 稳定边界恢复。
保留原目标函数、约束、接受规则；失败的 bootstrap 响应不丢弃、不重抽。
每次 exact 区间调用的内部时限仍为 10800 秒（3h），独立于作业时长。

资源与目录
按用户最新要求，每个计算任务申请：单节点、单任务、64 核、10G 总内存、03:00:00。
每个数组元素运行某一情形全部 1000 seed，不分块；数组没有 % 并发限制。
默认 1 个数组，8 个计算任务；另有 afterany 汇总任务，1 核、8G、1h，沿用原包设定。
三小时为申请上限，不保证所有情形一次完成；已处理 seed 可续跑。
线程配置：每个 seed worker 的 BLAS/OMP 为单线程，R 按 Slurm 分配核数并行 seed。
默认账户：def-liteep。
默认独立结果根目录：
${SCRATCH:-/scratch/$USER}/RD_v2_noneMAIC_MC1000_results_20261025
不影响原 full24 结果目录或已提交的 none+unit / score 作业。
环境变量使用独立前缀 RDNM，旧 RDFULL_NONE_* 设置不会覆盖本包资源。
可用 RDNM_ACCOUNT、RDNM_OUT、RDNM_CPUS、RDNM_MEM、RDNM_TIME 显式覆盖。

本机 PowerShell：上传与登录
$pkg = "D:\codexresults\ITC\rorqual\compareITC\anchored\RD_v2_noneMAIC_MC1000_20261025.tar.gz"
scp -o ConnectTimeout=60 -o ConnectionAttempts=3 "$pkg" "${pkg}.sha256" yanjin41@nibi.alliancecan.ca:~/
if ($LASTEXITCODE -ne 0) { throw "RD none+MAIC 上传失败" }
ssh yanjin41@nibi.alliancecan.ca

Nibi：校验、解压、提交全部 8 种情形
cd "$HOME" &&
sha256sum -c RD_v2_noneMAIC_MC1000_20261025.tar.gz.sha256 &&
tar -xzf RD_v2_noneMAIC_MC1000_20261025.tar.gz &&
cd "$HOME/RD_v2_noneMAIC_MC1000_20261025" &&
bash submit.sh

已上传旧 4h 包但尚未提交时，无需重新上传，可显式覆盖提交时限：
RDNM_TIME=03:00:00 bash "$HOME/RD_v2_noneMAIC_MC1000_20261025/submit.sh"
已提交的计算数组可在 Nibi 修改（第一个 ID 为计算数组，后面的汇总任务不改）：
out="${SCRATCH:-/scratch/$USER}/RD_v2_noneMAIC_MC1000_results_20261025"
job_id=$(head -n 1 "$out/active_job_ids.txt")
scontrol update JobId="$job_id" TimeLimit=03:00:00
squeue -j "$job_id" -o "%.24i %.20j %.10T %.6C %.10l %.8m %R"

无需从 Rorqual 迁移数据，也不依赖旧包目录。
脚本加载 StdEnv/2023 gcc/12.3 r/4.3.1；提交前只校验和准备配置，不做模拟。
查看本包作业：squeue -u "$USER" -o "%.24i %.22j %.10T %.6C %.10l %.8m %R"

结果与续跑
cells/<cell>/by_seed/rd_seed_00000001.rds 等是独立逐 seed 检查点。
每组保留 anchored_rd.rds、anchored_rd_estimates.csv、anchored_rd_coefficients.csv、
原方法 summary/status 表，以及 missing_seeds.csv、incomplete_seeds.csv。
report/progress.csv：各组记录数、已处理数、未完成数和缺失数。
report/all_methods_summary.csv：有效区间数、失败数、coverage、MCSE、宽度和 bias。
coverage_valid_pct 分母为有效区间数；returned_and_covered_pct 分母固定为计划 1000。
BC_boundary_truncated 是有效区间结果，单独计数并保留。
CASE_PROCESSED 表示 1000 条模拟记录处理完，不代表所有方法均返回有效区间。

若超时，待本包计算与汇总任务结束后，执行：
bash "$HOME/RD_v2_noneMAIC_MC1000_20261025/submit.sh"
会跳过处理完成的整组，并在其余组内复用同版本完成 seed。
未完成 seed 的旧记录先复制到 partial_history，再从其原始起点重新计算。
已完成但某个方法失败的记录保留，不自动删除、重抽或反复计算。
已有活动作业时拒绝重复提交；旧 none+unit 记录不能作为 none+MAIC 结果复用。
只续跑指定情形时，例如：bash submit.sh 1,3,8。
详细本地检查记录见 VERIFY.txt；不包含新的 coverage 试验或集群正式运行。
