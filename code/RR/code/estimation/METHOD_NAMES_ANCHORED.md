# Anchored 方法名称

2026-09-30 统一方法标签；anchored 表示间接比较类型，unit/maic 表示模型权重，二者独立。
默认 unit 模式中，旧 brm 与 brm_an 的点估计和方差相同。现仅保留单个 brm，
该方法遵循 model_weighting 配置；maic 模式也不再额外计算单位权重 BRM 基线。

主比较返回 18 个方法/推断标签：

`brm`, `LB`, `LP`, `RLP`, `CMH`, `CMH_sandwich`, `brmad`,
`brm-BC`, `brm-FC-BC`, `brm_b-BC`, `brm-FC`, `RLB`,
`LB-FC`, `RLB-FC`, `LP-FC`, `RLP-FC`, `LP-BFC`, `RLP-BFC`。

集群返回 18 条记录：brm、brm-FC、brmad 各 Wald/Exact 两条，其余 12 种各 Wald 一条。
含 EM 时 brmad / brm_b-BC 保留不适用状态；run_exact=FALSE 时保留 BC_not_run 状态。

| 旧标签 | 新标签 |
|---|---|
| `RLB_sandwich_an` | `RLB` |
| `CMH_sandwich_an` | `CMH_sandwich` |
| `RLB-FC_an` | `RLB-FC` |
| `RLP-FC_an` | `RLP-FC` |
| `RLP-BFC_an` | `RLP-BFC` |
| `LB-FC_an` | `LB-FC` |
| `LP-FC_an` | `LP-FC` |
| `LP-BFC_an` | `LP-BFC` |
| `brm_an` | `brm` |
| `brmad_an` | `brmad` |
| `LB_an` | `LB` |
| `LP_an` | `LP` |
| `RLP_an` | `RLP` |
| `CMH_an` | `CMH` |

RLB 使用原 RLB_sandwich_an 的拟合和方差，没有切换到已退役的旧 RLB_an 算法。
旧 brm 的独立 unit 基线已合并，brm-BC / brm-FC-BC / brm_b-BC 的名称保留。
新输出不使用 _an 或 _un 后缀，代码和主模拟输出文件名沿用 anchored 标识，
如 anchored_comparison.R、anchored_glm_fc.R、anchored_joint_estimates.csv。
variables/settings 另含 comparison_type="anchored" 和 method_naming_version=2。
新集群 checkpoint 的 settings.version=4，避免复用旧方法表。
历史数据不重写，独立绘图脚本同时识别新旧 CMH 名称。

## 本次数值核验

R 4.3.2 下，改名前后 unit/MAIC、含/不含 EM 四种设置中的系数与边际
点估计、SE、区间及状态在 1e-9 容差内一致。83 个 R 文件解析成功。
以下 13 项检查通过：check_unit_model_weights、check_glm_fc、
check_glm_fc_pipeline、check_glm_fc_cluster、check_joint_pipeline、
check_comparison_covariates、check_failure_handling、check_iteration_archive、
check_bc_score_reference、check_marginal_pipeline、check_custom_scores、
check_rare_boundaries、check_webster_pipeline。
cluster 检查用调用记录器核验 Exact 权重转发；marginal pipeline 包含真实小预算 BC。

check_marginal_op 在修改前和修改后均未通过其原有参数回收容差断言，
未改动估计算法来放宽该断言。check_robust_logbinomial 所调用的旧独立函数
依赖 sandwich 包，本机未安装，未完成该项；主比较 RLB 和 FC 不依赖该包。
