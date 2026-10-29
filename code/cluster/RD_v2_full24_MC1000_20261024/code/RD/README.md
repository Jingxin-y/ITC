# 独立 anchored RD 流程

此目录新增 RD 数据生成、估计、模拟、checkpoint、集群矩阵和汇总绘图。原 RR 文件及公共求解器不修改，通过上一级 `anchored_v1` 的相对路径加载。依赖与当前项目相同：R 4.3 系列和 MASS；不需要额外安装 PropCIs、sandwich、lmtest 或 logistf。

## 直接加法风险生成

借鉴 [Webster-Clark 等的 RD 模拟](https://academic.oup.com/aje/article/193/8/1176/7646053)，保留本项目三个 EM、三个预后非 EM，以及固定 BC 协变量概率、ESS 和精确 1:1 随机分配。这是项目改编，不是论文参数的原样复现。

```text
S_E = (1-z1)+z2+(1-z3)
S_P = z4+(1-z5)+z6
h_E = (S_E-.9)/sqrt(.63)
p_C = mu_C + beta_E*(S_E-.9) + beta_P*(S_P-.9)
p_A = p_C + (mu_A-mu_C) + lambda_AC*h_E
p_B = p_C + (mu_B-mu_C) + lambda_BC*h_E
```

生成不调用 OP/GOP 反解，不裁剪风险。当前为 **RD v2**：以 BC 人群为基准，将边际 A-B RD 设为 common 14 个百分点、rare 5 个百分点。C、B 保留原 RR 的 BC 均值；仅将 A 风险在全部协变量格上增加同一个常数。下面按百分比/百分点保留两位小数，内部继续使用校准精度。

| 场景 | C | A | B |
|---|---:|---:|---:|
| common | 41.86% | 64.22% | 50.22% |
| rare | 4.88% | 10.89% | 5.89% |

| 预设参数（百分点） | common | rare |
|---|---:|---:|
| beta_E | -16.85 | -2.00 |
| beta_P | 1.64 | 0.20 |
| lambda_AC low | 4.52 | 0.54 |
| lambda_AC high | 18.10 | 2.15 |

`lambda_BC=.5*lambda_AC`；high/low 比为 4。参数不依赖 ESS、来源、n 或 seed。`ess_source` 支持 `diffuse/em_driven/nonem_driven`，默认 ESS=.5。预设格概率满足 common `[.05,.95]`、rare `[.005,.995]`；自定义 `parameters` 仍须使所有三臂格风险严格位于 `(0,1)`，但不强制默认均值/包络，校准元数据会记录覆盖参数和是否匹配。

校准参考网格为两种事件率、两种强度、三种来源、ESS `.15/.25/.5/.8/1`。`data/rd_rr_calibration_comparison.csv` 包含全部 360 个逐人群逐臂对照值及 v2 有意增加的 A 风险；该增加量为 common 12.75、rare 4.81 个百分点。扣除这项增加后，原 AC-RR 最大差仍为 common 约 6.57、rare 0.77 个百分点。v2 不再要求 A 边际概率接近原 RR。

`data/rd_population_truth.csv` 包含 120 行逐场景逐人群的三臂风险及 A-C/B-C/A-B RD，`rd_population_truth_table()` 可重新生成。在默认 ESS=.5、diffuse 下：

| 场景 | AC A-B RD low | AC A-B RD high | BC A-B RD low/high |
|---|---:|---:|---:|
| common | 15.49 个百分点 | 19.97 个百分点 | 14.00 个百分点 |
| rare | 5.18 个百分点 | 5.71 个百分点 | 5.00 个百分点 |

BC 的 A-C RD 同时变为 common 22.36、rare 6.01 个百分点。其他 ESS/来源下 AC 的 A-B RD 仍随 EM 人群分布变化，common 在强偏移场景可超过 20 个百分点；10%-20% 是本次确认的默认场景范围。

## 方法、尺度和不确定性

主比较目标是 **BC 人群中的边际 A-C RD**，与原 RR 主流程的人群/治疗目标一致。本次不新增 A-B 间接比较的估计和 CI；B-C、A-B 总体真值另存。

九条方法/推断记录：`brm`、`brm-BC`、`brm-FC`、`brm-FC-BC`、`IB`、`RIB`、`RLPM`、`MN_eff`、`RD_MAIC_sandwich`。不加入 Bayes/GC，也不套用原 log-link GLM FC/BFC。

`estimation_mode="scores"` 是原默认模型；`"none"` 不在结局模型中纳入协变量：BRM/FC 的 RD 与 OP 都只有截距，IB/RIB/RLPM 使用 `y ~ x`。MAIC 校准变量由 `weight_covariates` 独立控制，`all` 仍校准原始六个协变量；MN_eff 与加权臂均值 sandwich 不受估计模式影响。生成数据及总体 RD 真值不变。结果、诊断和 checkpoint 保存估计模式并拒绝混用。

- BRM/FC 复用公共 RD 求解器。当前默认 `brm_effect_model="linear_atanh"`：`atanh(RD)=alpha_0+alpha_1*h_E`，效应部分两个参数；OP 仍按所选 EM 总分的取值饱和。`h_E` 使用 BC 样本均值及标准差构造，生成数据不变。此模型属于较简约的效应模型，通常不能精确表达当前自然 RD 线性真值，尤其 common/high 下需要评估模型偏离。
- 显式 `brm_effect_model="saturated"` 可复现原来的四层饱和效应模型（截距加三个 dummy）；它能准确表示当前加法 RD 真值。内部 alpha 都是 atanh(RD) 的基底系数，不能直接当作边际 RD；两种模型都计算 MAIC 标准化后的平均 `tanh(va %*% alpha)`。`estimation_mode="none"` 继续使用截距模型。
- RD 的 BRM/FC 点估计、导数、固定/联合方差、边际风险及 BC 反演统一使用 `rd_probabilities()` 的稳定 RD/OP 转换。`rd_solver_adapter.R` 在私有环境中复用共享求解器，内部 `getProbRD/getProbScalarRD` 绑定到稳定转换，不覆盖全局 RR/公共函数。移除旧 OP +/-12 的 0/1 近似，并处理零 RD、极小 OP 的下溢。共享求解器原有内部概率界限、收敛检查及失败状态保留；极端有限预测值仍可能因机器精度舍入到 0/1。
- IB 使用 binomial identity link；RIB 共用该拟合，使用 observed-bread sandwich。RLPM 使用线性概率模型及 HC3。默认回归包含 `h_E/h_P` 主效应和 `x:h_E`；先生成两组反事实预测，再用 MAIC 权重求平均差。
- 默认 `model_weighting="unit"`：拟合权重全 1，标准化权重仍为 MAIC，方差固定标准化权重与协变量。`"maic"` 对 BRM/FC、RIB、RLPM 和加权 RD 计算联合方差，包含权重、经验标准化及 BC 协变量均值不确定性；IB 保留固定权重 model-based 方差。FC 沿用现有 RR 的一阶未修正 score 方差约定。
- MN_eff 每臂采用 Kish `n_eff=(sum(w))^2/sum(w^2)` 与加权事件率进行 MN score 反演，保留 MN 总有效人数的 `N/(N-1)` 修正。这是固定 MAIC 权重的有效计数近似，不是标准整数计数 MN 的精确延伸。[标准 MN 文档](https://search.r-project.org/CRAN/refmans/PropCIs/html/diffscoreci.html)。`se=NA`，区间宽度诊断另放 `se_from_ci_width`，不计入 mean model SE。
- RD_MAIC_sandwich 与 MN_eff 共用每臂加权事件率的点估计，直接用加权臂均值的 influence function 计算方差；不使用分层共同 RD 假设。
- RD BC 沿用当前 RR 的候选值与固定观察目标双重 bootstrap/Blaker 型反演，区间与 P 值均在自然 RD 上。给定协变量、治疗分配、模型和标准化权重；不会重新抽样 MAIC/BC 汇总，不保证有限样本精确覆盖率。`run_exact=FALSE` 保留 `BC_not_run` 状态。

Wald 区间不裁剪至 `[-1,1]`，LPM 的反事实风险不裁剪。`counterfactual_risks_valid` 与点估计、SE、CI/P 状态分别保存；单个方法失败不会删除其他方法或独立可计算的点估计。未出现/单臂缺失的分数格可能使饱和 BRM 不可识别，其状态保留。

## 本地使用

在 `anchored_v1` 目录运行：

```r
source("RD/load_anchored_rd.R", encoding="UTF-8")
source_anchored_rd()
design <- webster_rd_scenario_design("common", .5, "diffuse", "high")
result <- run.anchored.rd(design, 500, seed=19, run_exact=FALSE)
result$estimates
result$coefficient_estimates
design$marginal_truth
```

已有数据可调用 `get.compare.anchored.rd(data.ipd,data.ad,...)`。可设置 `score_spec`、`weight_covariates` 和回归的 `model_covariates/effect_modifiers`；`brm_effect_model` 独立控制 BRM/FC 的效应基底，不改变其他回归模型或 MAIC 校准集合。`data.ipd` 是包含 `data` 的 list，`data.ad` 使用原 `get.aggre()` 汇总格式。

直接拟合 RD BRM/FC 时使用 `rd_get_estimate(..., param="RD")`；直接转换概率使用 `rd_getProbRD(theta, phi)`（返回矩阵）或 `rd_getProbScalarRD(theta, phi)`（返回长度二向量）。这些入口支持标量/向量、两列预测矩阵和标量广播。公共 `get.estimate()/getProbRD()` 继续保留原实现，RD 主流程通过专用适配调用稳定版本。`result$variables$probability_mapping` 及 BRM/FC fit 的同名字段标记为 `stable_rd_op_v1`。

```sh
Rscript RD/estimation/run_simulation_rd.R RD/simulation_config_rd.R RD/results/local
Rscript RD/summary/summarize_rd.R RD/results/local/anchored_rd.rds
Rscript RD/summary/plot_rd_results.R RD/results/local
Rscript RD/tests/check_rd_pipeline.R
Rscript RD/tests/check_rd_failure_handling.R
Rscript RD/tests/check_rd_stable_solver.R
Rscript RD/tests/check_rd_estimation_modes.R
Rscript RD/tests/check_rd_effect_models.R
```

示例配置为 common/high/diffuse、n_AC=n_BC=1000、10 个 seed、unit 拟合、BC 关闭。输出 `anchored_rd.rds`、`anchored_rd_estimates.csv`、完整系数 CSV、汇总/方法状态/参数状态 CSV。checkpoint 要求相同 RD 设计、分析参数、seed 序列及代码哈希，否则拒绝复用；v1 checkpoint 不能与当前 v2 合并。CLI 最后一个可选参数为公共代码目录。

稳定求解器接入前的 v2 checkpoint 也不能复用；数据生成仍是同一 v2，概率转换标识及代码哈希会区分分析实现。旧验证结果保留在 `validation/v2_before_stable`。

汇总复用原尺度无关函数，RD adapter 增加 RMSE、区间失败率和包含失败重复的覆盖率。图中覆盖率采用有效 CI 条件下的覆盖率，失败比例独立报告；正式研究应同时查看 CSV 中的 `coverage_all_replicates` 和 `interval_failure_rate`。

## 集群入口

`submit_rd_diffuse_mc1000.sh` 对应 RR 的 `submit_diffuse_mc1000.sh` 分组。all/nonem 每组 24 个场景（common/rare、low/high、n=100/500、ESS=.15/.5/.8）；em 每组 48 个场景（n=100/200/500、ESS=.15/.25/.5/.8）。每个场景 1000 次外层重复、seed=1:1000、full CI。scores 使用 unit 拟合、32 核、12 小时；none 使用 MAIC 拟合、16 核、3 小时；两组均为 15 GB。

```sh
export RD_CODE_DIR="$HOME/MAIC_ess_em_20261004"
export RD_RESULT_ROOT="$SCRATCH/MAIC_RD_MC1000_20261007"
bash "$RD_CODE_DIR/RD/submit_rd_diffuse_mc1000.sh" all scores dry-run
bash "$RD_CODE_DIR/RD/submit_rd_diffuse_mc1000.sh" all none dry-run
bash "$RD_CODE_DIR/RD/submit_rd_diffuse_mc1000.sh" all scores submit
bash "$RD_CODE_DIR/RD/submit_rd_diffuse_mc1000.sh" all none submit
```

该入口与更新后的 `submit_rd_grid.sh` 配套，后者支持 `RD_RESULT_ROOT`；未设置时单场景仍使用 `$SCRATCH/MAIC_RD`。部署时用 `RD_MC1000_UPDATE_SHA256SUMS` 校验两个 shell 文件及此前 none 更新的三个 R 文件。MC1000 指外层模拟次数，full BC 的 bootstrap 预算沿用现有值。

```sh
bash RD/submit_rd_smoke.sh                 # dry-run，不提交
bash RD/submit_rd_matrix.sh               # 144 场景 dry-run
bash RD/submit_rd_smoke.sh submit         # 明确提交小规模 smoke
bash RD/submit_rd_grid.sh submit 200 0.5   # 明确提交单个 full 场景
bash RD/submit_rd_matrix.sh submit        # 明确提交 144 场景
bash RD/submit_rd_grid.sh collect /path/to/one/run
```

支持 `RD_EVENT/RD_ESS_SOURCE/RD_EM_LEVEL/RD_WEIGHT_MODE/RD_MODEL_WEIGHTING`，`RD_CI_MODE=full/smoke/wald`、`RD_EST_MODE=scores/none`。`RD_BRM_EFFECT_MODEL=linear_atanh/saturated` 控制 BRM/FC 效应模型，默认 `linear_atanh`，结果文件夹含 `brm_linear_atanh`；none 模式标记为 `brm_intercept_only`。默认 full 为 500 重复、32 cores、10 GB、3 小时、账户 `def-liteep`；smoke 为 2 重复/2 cores/30 分钟。可用 `RD_NSIM/RD_CORES/RD_FIRST_SEED/RD_ACCOUNT/RD_MEM_GB/RD_WALLTIME/RD_R_MODULE` 覆盖。smoke 的 bootstrap 预算只检验运行路径，不能用于正式覆盖率结论。

例如 `RD_EST_MODE=none bash RD/submit_rd_grid.sh submit 200 0.15`，场景名含 `est_none`。旧入口不识别此变量；需同时更新 `submit_rd_grid.sh`、`estimation/rd_cluster.R`、`estimation/rd_simulation.R` 与 `estimation/rd_comparison.R`。已经提交的任务使用自己的代码 snapshot，后续上传/环境变量变更不会改变其资源或模型。

每个 seed 保存独立 RDS，采集器拒绝混合 RR、不同设计/分析或源代码版本。结果存于 `$SCRATCH/MAIC_RD`。提交时保存完整公共代码及 RD 代码 snapshot；正式任务未随实现自动提交。

线性 h_E 更新需要同时上传 `estimation/rd_comparison.R`、`estimation/rd_simulation.R`、`estimation/rd_cluster.R` 和 `submit_rd_grid.sh`。绘图脚本 `summary/plot_rd_cluster_results.R` 可同时更新，以在图和 metrics 中标记效应模型；旧结果缺少该设置时按原饱和模型显示。检查点记录模型设置及源哈希，旧版本不能直接续算，linear_atanh 与 saturated 不能混合采集。此次不改变 RD/OP 概率映射、FC 导数、BC 算法和预算、数据生成、权重求解、其他方法或 RR/公共文件。

## 验证记录

`VALIDATION.md` 记录实际执行的验证与限制。测试使用 base R/MASS，输出仅写入临时目录或 `RD_TEST_OUTPUT`；不改动代码或原 RR 结果。
