## 在本版本工作目录逐段运行：2026-10-15 RR effect v2。
source("load_anchored.R", encoding = "UTF-8")
source_anchored_data()

## 1. 修改这一行即可切换场景；Common/common/rare/Rare 均可。
## v2 固定 BC 均值：common C/A/B=.418596/.514669/.60；rare=.048831/.10/.15。
event <- "Common"  # 改为 "rare" 使用罕见事件场景
design <- webster_rr_scenario_design(
  event = event,
  ess_ratio = 0.50, ess_source = "diffuse", em_level = "high"
)
## common/rare 分别校准 b、治疗及预后乘数，baseline_risk 仍不是总体平均风险。
## 默认两项治疗不同；改变某个乘数不再保证原事件率，且会检查概率是否合法。
## 改变 ESS/q 时风险表不变，AC 事件率可随人群分布改变。
design$event               # 规范化后的场景名
design$parameters          # 直接风险乘数，不是 OP/GOP 参数
design$covariate_roles     # z1:z3 是 EM，z4:z6 是非 EM
design$coefficient_truth   # 固定总体 h_E 尺度下的条件 alpha/lambda
design$e                   # 固定的 AD / BC 概率
design$q                   # 根据 ESS 与 direction 求出的 IPD / AC 概率
design$ess_ratio           # 理论相对 ESS，不是某次样本实现的 ESS
design$risk_table          # 两试验共用的 16 行风险表
design$marginal_truth      # 分别在 AC、BC 人群下积分的风险及边际 log RR
design$calibration         # v2 BC 事件率目标及校准结果

## 可独立检查任意分数组合的三臂风险，无需 OP/GOP 求根。
tab <- design$risk_table
cbind(p_C=tab$p_C, p_A=tab$p_A, p_B=tab$p_B)

## 以下 webster_rr_design() 是保留的旧参考设计，不使用上面的 v2 目标。
## 旧参考设计可直接输入任意六维 q，此时不按 ess_ratio 反解：
# design <- webster_rr_design(
#   event=event, q=c(.5,.5,.5,.5,.5,.5)
# )
## 可显式覆盖单个系数；其他系数继续使用所选场景值：
# design <- webster_rr_design(event=event, rr_BC=.80, em_rr_BC=.90, ess_ratio=.50)
## 保留原 RR-GOP 对照时显式使用：
# old_design <- joint_rr_gop_design(event="common", ess_ratio=.50)

## 2. 同一 design 生成两项独立随机试验；保持项目原来的 z1,...,z6 定义。
## n_AC、n_BC 均为正偶数；本例每项试验的两组各有 500 人。
set.seed(1)
sim <- generate_anchored_data(n_AC = 1000, n_BC = 1000, design = design)
data.IPD <- sim$data.IPD
data.AD.or <- sim$data.AD.or
data.AD <- sim$data.AD
head(data.IPD$data)
data.IPD$count
data.AD.or$count

## 3. 查看 AD 汇总；get.aggre() 在 data/ad_summary.R 中定义。
data.AD$anchored.mean.x
data.AD$anchored.cov.x
data.AD$anchored.n
## 若后续按分数列求权重，先汇总对应的列：
# score.AD <- get.aggre(data.AD.or, covariates = "h_E")

## 求权重的示例见 estimation/example_weights.R。
## 旧版比较见 estimation/legacy/example_anchored_legacy.R；
## 旧 alpha 真值不可用于本例的边际 RR 分析。
