# RD v2 与稳定求解器验证（2026-10-09）

## 对应 RR 分组的 RD MC1000 入口

新增 `submit_rd_diffuse_mc1000.sh`，沿用 RR 分组的样本量、ESS、common/rare、low/high、1000 次外层重复及 full CI：scores/unit 为 32 核/12 小时，none/maic 为 16 核/3 小时，内存按用户要求统一 15 GB。`submit_rd_grid.sh` 仅增加 `RD_RESULT_ROOT`，未指定时保留原默认路径；公共/RR 源代码和 RD 估计 R 文件未改动。

两组 all 的 24 场景 dry-run 校验通过；本地假 sbatch 验证实际提交参数的核数、内存、时限和结果路径。实际 none/MAIC 的单 seed 小预算 RD smoke 返回九条可用点估计、CI/P 值，覆盖联合方差与真实 BC 反演。测试输出位于 `validation/mc1000_update`；没有向集群提交正式作业。

## all 校准 / none 估计与资源更新

新增 `estimation_mode="scores"/"none"` 和 `RD_EST_MODE`。none 使 BRM/FC 的效应与 OP 仅含截距、IB/RIB/RLPM 使用 `y ~ x`；all 仍校准六个原始协变量。默认 full 资源按用户要求改为 32 核、10 GB，重复数 500 与时间 3 小时沿用原值。模型、并行 worker、结果列、诊断、场景名和 checkpoint 一致记录估计模式。

实际检查：生成原始数据/BC 原始汇总/MAIC 权重在 scores 和 none 下完全相同；none 设计矩阵与独立臂均值对照正确；MN/加权臂均值方法不随估计模式变化；非法模式拒绝；none checkpoint 可恢复且拒绝混入 scores；两个 worker 实际执行 none；真实小预算 none BC smoke 返回九条可用 CI/P。原默认 scores 的点估计、方差、状态和系数与更新前逐项完全相同。用户指定的 common/rare、n=200、ESS=.15/.8、high/diffuse 四个 none 场景在 seed=19 均返回九条可用点估计。

完整 RD pipeline 回归检查通过，包括原 scores 模型、联合方差、真实 BC、检查点、汇总/绘图及两个 worker。八个 scores/none 提交 dry-run 均显示 32 cores、10G，Bash 语法和非法模式拒绝通过。验证记录位于 `validation/none_update`，未提交正式集群作业。早期验证结果为当时源版本，源代码哈希会阻止直接复用其 checkpoint。

`tests/check_rd_estimation_modes.R` 的 BRM 独立臂均值对照显式使用严格优化容差；生产默认容差及共享求解器保持原值。原共享/RR 文件未修改。

## 稳定概率转换接入

RD 主流程的 BRM/FC 点估计、FC 导数、固定及联合方差均通过 `rd_solver_adapter.R` 的私有环境接入稳定 RD/OP 转换，与边际标准化和 RD BC 反演一致。全局 RR/公共函数不被覆盖，原文件不修改。诊断标识为 `probability_mapping="stable_rd_op_v1"`。

本次实际完成：

- 新增 `check_rd_stable_solver.R`：全部标量/向量/两列矩阵/成对输入/广播/空输入接口；非法输入拒绝；零 OP、接近零 OP、旧 +/-12 阈值附近及极端有限预测值。
- 映射与独立风险空间 root oracle 一致；满足 RD 及可分辨内部风险的 log-OP 约束。
- 移除旧 OP 阈值的 0/1 近似；处理 RD=0、log-OP=-1000 的 0/0，以及极小非零 RD 下指数下溢但风险仍可表示的情况。
- FC 一阶、二阶导数与稳定映射的独立数值差分/隐式约束求导一致。
- 临时使全局旧 `getProbRD/getProbScalarRD/firth_rd_derivatives` 报错，实际 unit/MAIC 两条 BRM/FC 流程仍返回九个可用点估计，证明新路径已经调用。
- 加载 RD 适配前后，全局 RR/公共函数体、参数及环境保持一致，RR 概率结果一致。
- 完整 RD pipeline 和失败处理检查通过，包括两个 worker、检查点恢复、汇总/绘图和真实小预算 BC 反演。
- 当前稳定版本 n=500、common/high/diffuse、seed=19 的集群式 smoke 返回九个可用区间/P 值，记录在 `validation/cluster_smoke`。
- 生成数据和真值与接入前完全相同；实际旧求解版本的 v2 checkpoint 被新版本拒绝。
- 在加载 RD 适配后重新运行五项原 RR/公共回归检查，全部通过；127 个原文件的 SHA256 全部保持一致。

当前验证记录在 `validation/check_rd_stable_solver.txt` 与 `validation/stable_integration_checks.txt`。接入前的 v2 结果与说明保留在 `validation/v2_before_stable`。共享求解器原有概率内部界限、方差保护及失败状态保留；极端有限预测值可能因机器精度舍入到 0/1。没有执行正式大规模模拟或 SLURM 提交。

## 数据生成 v2

用户确认的目标：BC population 的边际 A-B RD 为 common 14.00 个百分点、rare 5.00 个百分点；默认 ESS=.5、diffuse 的 AC low/high 分别约为 common 15.49/19.97、rare 5.18/5.71 个百分点。

只将 A 的每个协变量格风险增加 common 12.75、rare 4.81 个百分点。C、B 风险、协变量分布、ESS 来源、high/low 四倍比和 B/C 修饰强度比 .5 均保持原值。BC 的 A-C RD 随之为 common 22.36、rare 6.01 个百分点。

## 本次实际通过的检查

- 60 个生成设计、每个设计的全部 64 个协变量组合：风险合法，common 在 [.05,.95]、rare 在 [.005,.995]。
- 逐格核对只有 A 增加同一个常数；C、B 均未改变。
- 120 行人群真值：所有场景的 BC A-B RD 精确为 .14/.05，rare 的 AC A-B RD 不低于 .05。
- 默认 diffuse/ESS=.5 的 common A-B RD 在 [.10,.20]。完整参考网格的 AC A-B RD 范围为 common 14.00-28.40、rare 5.00-6.71 个百分点；common 的范围要求适用于用户确认的默认场景。
- 360 行 RR 对照：BC 的 C/B 均值仍精确匹配原 RR；A 的差距显式记录为有意调整。扣除 A 的调整量后，原 AC-RR 校准界限不变。
- 完整 `check_rd_pipeline.R` 通过：可复现、饱和 BRM 表达、数值梯度、MN 独立 oracle、HC3/identity 模型/加权 RD sandwich、unit/MAIC 联合方差、真实小预算 BC 反演、checkpoint 恢复、汇总绘图和两个 worker。
- `check_rd_failure_handling.R` 通过：零事件、缺失 EM 格、错误权重/尺度/BC 参数、BC timeout、profile 梯度及 rare/低 ESS 的独立状态。
- 将实际 v1 单 seed checkpoint 交给 v2 流程，按设置不兼容拒绝复用。
- v2 实际集群式 smoke：common/high/diffuse、n=500、seed=19，九条记录的 CI/P 状态均为 available，结果位于 `validation/cluster_smoke`。
- 127 个原 RR/公共文件 SHA256 与实施前记录完全一致。原回归检查和 144 场景 dry-run 的 v1 执行记录保留在 `validation/v1/VALIDATION.md`；本次不改动其代码及提交脚本。

## 文件与限制

`data/rd_population_truth.csv` 是当前 v2 总体真值；`data/rd_rr_calibration_comparison.csv` 是当前 RR 对照。`validation/v1` 保留旧版本的结果和说明，旧结果不代表当前生成过程，不能与 v2 混合。检查点中的 DGM 版本和源代码哈希会阻止混用。

执行环境为 Windows/R 4.3.2/MASS 7.3-60。Windows 入口使用可用的 UTF-8 locale；本地并行 socket 在沙箱外验证。正式 SLURM 任务未提交，大规模 Monte Carlo 未执行。smoke 的极小 bootstrap 预算仅验证运行路径，不能用于正式覆盖率结论；两重复绘图同样仅验证输出。
## 2026-10-13：线性 h_E 的 atanh-RD 效应模型

按用户要求，scores 的 BRM/FC 默认改为 `atanh(RD)=alpha_0+alpha_1*h_E`；OP 仍按 EM 总分取值饱和。显式 `brm_effect_model="saturated"` 保留旧模型，none 仍为截距。RR/公共代码、稳定概率映射、FC 导数、BC 求解算法/预算、数据生成和 MAIC 权重算法保持原值。

实际通过的验证：

- 默认两列效应矩阵、四列 OP；旧饱和四列效应和 none 截距均可调用。数据、BC 汇总、MAIC 和其他五个方法不随效应模型变化。
- 线性模型的边际 RD 和全参数梯度与独立加权平均/数值差分一致。
- 之前 5 个缺臂 seed（11、30、33、37、58）在 rare/ESS=.15/diffuse/n100 的 low/high 下，两种线性 BRM/FC 全部返回点估计；显式饱和模型仍保留原不可识别失败。
- 本地用 seed=1:100 重跑 low/high 的点估计：每个情形 BRM 与 FC 均返回 89/100；11 次仍为 weighting_failed，没有新增拟合或目标失败。此次 100 seed 检查关闭 BC，不代表区间返回率或覆盖率。
- common/rare、scores/none、n500/seed19 的显式 saturated 结果与更新前的点估计、SE、状态和系数逐项相同；none 的默认结果也保持相同。
- 新模型 checkpoint 恢复通过，切换模型拒绝；采集器拒绝混入 saturated，两个本机 worker 正确传递并运行 linear_atanh。
- 一个 previously-unidentified rare n100/ESS=.15/seed11 的小预算 BC 实际完成，BRM/FC 都返回自然 RD 区间，使用原多初值/nlminb rescue。完整 RD pipeline、稳定概率/FC 导数、方法失败处理及原 RR 重试检查通过。
- 144 场景矩阵和 MC1000 scores 的 24 场景 dry-run 通过，资源设置沿用现值；未提交正式集群作业。更新后的独立绘图脚本正确标记模型并输出 PNG/PDF。
- 核对 111 个公共源文件的 SHA256，全部未改变。BC 源文件本身未修改。

验证脚本及输出位于 `anchored/RD_linear_hE_update_20261013/validation`。本机受限环境的并行 worker 连接检查在允许本机套接字的执行模式下完成；第一次完整 pipeline 权限审核超时，重试后通过。

限制：自然 RD 仍按当前加法公式生成，atanh-RD 线性模型通常存在模型偏离，尤其 common/high。尚未用正式 full BC 或大规模模拟验证新模型的偏倚、覆盖率或区间返回率。MAIC 不可行与搜索截断并没有通过此次模型修改得到解决。

