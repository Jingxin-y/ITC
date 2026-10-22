RR effect v2 — 2026-10-15

本版本以 anchored_v1 在本次修改开始时的当前 RR 代码为基础，包含 data、
estimation、summary 及 ESS 集群提交入口。SOURCE_MANIFEST.csv 保存修改前
所复制文件的 SHA256。原目录和已有模拟结果未被改写；RD 子项目不在此包中。

1. BC 目标总体的平均风险（C / A / B）

common: 0.418596226792393 / 0.514669087904676 / 0.600000000000000
rare:   0.0488309131792698 / 0.100000000000000 / 0.150000000000000

rare 的对应关系明确为 A=.10、B=.15。common 的 A 和两种事件的 C 保留原值。
以上都是 BC 总体真值，不是每个样本必须达到的事件率。AC 总体平均风险继续
随来源分布变化。low/high EM 和各 ESS 来源/大小下，BC 的这些目标保持一致。

以 A/B 为方向：common RR_AB≈.85778、RD_AB≈-.08533；
rare RR_AB=2/3、RD_AB=-.05。实际研究仍需区分估计的 A/C 目标与完整 A/B ITC。

2. 为保证概率合法而调整的参数

原 baseline_rr_E=.55 在 common、low EM、B=.60 时会使某个评分格风险约为
1.028，无法生成合法 Bernoulli 结局。现统一改为 .59；baseline_rr_P 保持 1.08。
全部 16 个评分格风险严格检查在 (0,1)，不裁剪风险。该改变也会改变条件风险
形状与 RD 异质性，因此本版本属于新的 DGM，不应与旧结果混合。

EM 斜率保留：low lambda_AC=.15、lambda_BC=.075；
high lambda_AC=.60、lambda_BC=.30。三个 EM、三个非 EM、评分方向、
BC 协变量分布、ESS 定义和求解、1:1 分配及所有估计/CI 算法沿用复制时的代码。

新目标单独放在 data/generate_webster_RR_scenarios.R 的
webster_rr_scenario_targets()，不修改旧 webster_event_targets()。
webster_rr_design()、joint_rr_gop_design()、webster_rd_design() 仍是历史对照。
获得新目标必须调用 webster_rr_scenario_design()，或使用下面的新版本入口。

3. 在本目录运行 R

source("load_anchored.R", encoding="UTF-8")
source_anchored()
design <- webster_rr_scenario_design(
  event="rare", ess_ratio=.15, ess_source="em_driven", em_level="high")
design$marginal_truth
result <- run.anchored(design, n_AC=200, n_BC=200, seed=1,
  score_spec=make_brm_score_spec(), weight_covariates=paste0("z",1:6),
  model_weighting="unit", return_data=TRUE, run_exact=FALSE)

批量配置 estimation/simulation_config.R 默认采用 common/high/diffuse 的 v2 设计：
Rscript estimation/run_simulation_anchored.R estimation/simulation_config.R results_v2

Windows 如果 R 因不支持 C.UTF-8 而无法读取中文注释，可先在 R 会话执行：
Sys.setlocale("LC_CTYPE", "English_United States.utf8")

4. 集群

把此文件夹整体上传到一个新的代码目录；进入该目录后运行。
如当前 shell 曾设置 MAIC_CODE_DIR，须先更新为本目录，避免误用旧源码：

export MAIC_CODE_DIR="$PWD"
bash submit_ess_em_smoke.sh

正式单场景示例（提交，不是 dry run）：
MAIC_EVENT=common MAIC_EM_LEVEL=high MAIC_ESS_SOURCE=diffuse \
  MAIC_WEIGHT_MODE=all MAIC_MODEL_WEIGHTING=unit \
  bash submit_ess_em_grid.sh 200 0.5

先列出原 144 个实验单元：
bash submit_ess_em_matrix.sh dry-run

全部提交：
bash submit_ess_em_matrix.sh submit

新提交输出位于 $SCRATCH/MAIC_ess_em_v2/。默认模拟次数、核数、时间和
MAIC_MODEL_WEIGHTING 等参数沿用原提交脚本，可按原有环境变量设置。
主估计入口仍默认 unit 模型拟合；model_weighting="maic" 可显式选择。
CMH、数值修正和方差口径均来自本次复制的当前代码，不回退到较早快照。

版本标记：dgm=webster2024_rr_ess_em_v2；scenario version=2；
risk_target_version=BC_CAB_20261015；checkpoint settings.version=11。
不要把旧版本 seed 文件移入新结果目录；版本及 DGM 参数不同会拒绝复用。

5. 检查

Rscript data/tests/check_webster_rr_scenarios.R
Rscript estimation/tests/check_webster_rr_scenario_pipeline.R
Rscript estimation/tests/check_ess_em_cluster_settings.R
Rscript estimation/tests/run_cluster_ess_em_smoke.R . smoke_v2 1

前者以独立 64 个原始协变量组合积分检查 60 个 event/EM/source/ESS 设计，
核对新目标、合法风险、既有斜率与 ESS 分解。smoke 使用小 bootstrap 次数
验证完整执行链，不构成覆盖率验证，也不是正式大规模模拟结果。
