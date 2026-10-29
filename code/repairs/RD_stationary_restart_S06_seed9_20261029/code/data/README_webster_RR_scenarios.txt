ESS 来源和 EM 强度的 Webster RR 新 DGP（2026-10-04）

Material Passport
依据：用户已确定的 ESS 来源/EM 强度方案；固定当前 BC 平均风险。
范围：新 DGP、抽样入口、加载入口、使用示例和检查脚本。
该方案是本项目的扩展，并非论文情形的逐项复现。

加载：
source("load_anchored.R", encoding="UTF-8")
source_anchored_data()

创建总体设计：
design <- webster_rr_scenario_design(
  event="common", ess_ratio=.15, ess_source="em_driven", em_level="high")

抽样：
set.seed(1)
sim <- generate_anchored_data(200, 200, design)

也可以一次完成：
sim <- generate_webster_rr_scenario(
  n_AC=200, n_BC=200, event="common", ess_ratio=.15,
  ess_source="em_driven", em_level="high", seed=1)

n_AC/n_BC 是每项试验总人数，两臂各 n/2。必须是正偶数。
不提供 seed 时使用调用者的当前随机数状态。
返回 data.IPD、data.AD.or、data.AD、design，与原抽样器接口相同。

参数：
event：common / rare。
ess_ratio：Full 校准的总体理论相对 ESS，范围 (0,1]；包括 0.25。
ess_source：diffuse / em_driven / nonem_driven。
em_level：low / high。

来源定义：
diffuse：a=.5，r_E=sqrt(r)，r_P=sqrt(r)。
em_driven：a=1，r_E=r，r_P=1；只有 z1:z3 的 AC/BC 分布不同。
nonem_driven：a=0，r_E=1，r_P=r；只有 z4:z6 的 AC/BC 分布不同。
design$relative_ess 给出 full/em/non_em 的理论相对 ESS。
总体相对 ESS=1 不表示有限样本 MAIC 拟合权重必定全为 1。
比较校准集合时必须使用同一份数据，不能为了重新满足子集 ESS 而改变 DGP。

EM 参数以固定 BC 标准差 sqrt(.63) 为单位：
low：lambda_AC=.15、lambda_BC=.075、lambda_AB=.075。
high：lambda_AC=.60、lambda_BC=.30、lambda_AB=.30。
BC/AC 的系数比 rho=.5；原始评分乘数为 exp(lambda/sqrt(.63))。
这里的 low/high 是围绕原约 .30 系数的一半/两倍选定的实验水平，没有通用强弱阈值含义。

S_E=(1-z1)+z2+(1-z3)，S_P=z4+(1-z5)+z6。
pC=b*.55^S_E*1.08^S_P。
pA=pC*rr_AC*em_rr_AC^S_E。
pB=pC*rr_BC*em_rr_BC^S_E。
b、rr_AC、rr_BC 按 event/em_level 通过独立 Bernoulli 的解析乘积矩求解。
基线形状 .55/1.08 在所有新情形中一致，较原预设有所改变，以容纳较强 EM。

固定 BC 总体平均风险（C/A/B）：
common：.418596226792393 / .514669087904676 / .502197586315061。
rare：.0488309131792698 / .0608090036668075 / .0589376565366341。
AC 风险随来源分布改变，不额外匹配旧 AC 均值。
参数不随 ESS/source/n/seed 重新校准。所有 16 种评分组合的风险必须严格处于 (0,1)，无裁剪。
design$calibration 只描述 BC 均值匹配；design$scenario 保存版本、来源、强度和固定尺度。

与估计流程衔接：
source_anchored()
scores <- make_brm_score_spec(paste0("z",1:3), paste0("z",4:6))
result <- run.anchored(design, 200, 200, seed=1, score_spec=scores,
  weight_covariates=paste0("z",1:6), model_weighting="unit",
  bc_target="marginal", run_exact=TRUE)

模型拟合用全 1 权重，marginal 风险及对应 CI 目标用 MAIC 权重，由现有估计入口负责。
weight_covariates 可指定 Full=z1:z6、E-only=z1:z3、P-only=z4:z6。
完整 cluster 入口为根目录的 submit_ess_em_grid.sh；小测试为 submit_ess_em_smoke.sh。
旧 submit_common_grid.sh 等入口仍默认调用原 webster_rr_design。
R 配置脚本可以直接将 design 改为 webster_rr_scenario_design。

检查命令（在 anchored_v1 下）：
Rscript data/tests/check_webster_rr_scenarios.R
Rscript estimation/tests/check_webster_rr_scenario_pipeline.R
Rscript data/tests/check_webster_rr.R

前者检查 60 个 source/ESS/event/EM 组合、64 格总体积分、抽样接口和可复现性。
后者检查新 DGP 接入现有估计流程及权重设置，run_exact=FALSE。
这些是实现检查，不是正式覆盖率模拟，也不保证稀有事件样本总能返回点估计和 CI。
