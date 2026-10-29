# Webster-Clark RR 生成机制：三个 EM、三个非 EM

本实现借鉴 Webster-Clark et al. (2024) 的 **RR 敏感性分析**中直接建立
log-risk 模型的方式，保留本项目的三 EM / 三非 EM 结构。
它是 adapted DGM，不是对论文模拟的逐项复现。

来源：[论文](https://doi.org/10.1093/aje/kwae048)；
[作者 RR SAS 代码](https://github.com/mawcpharmdphd/ext_var_sel/blob/6d495d97b5002d16e5ae50bfd2f33c1a58f3dddc/Sensitivity_analyses/simulating_variable_selection_RR.sas)。
核对的仓库版本：`6d495d97b5002d16e5ae50bfd2f33c1a58f3dddc`。

论文原始 RR 模型为
`p(Y=1|X,Z)=.90*.67^(X+Z2+Z3+Z5+Z6+X*Z3+X*Z6)`。
原文的 EM 是 Z3、Z6，另有两个预后非 EM 和两个无关变量。
这里按用户要求保留 z1–z3 三个 EM、z4–z6 三个预后非 EM，
不采用原文的六类变量分组，也不采用额外含相关 Z7 的敏感性场景。

## 当前生成公式

```text
S_E = (1-z1) + z2 + (1-z3)
S_P = z4 + (1-z5) + z6
h_E = (S_E-.9)/sqrt(.63)
h_P = (S_P-.9)/sqrt(.63)

log p_C = log(baseline_risk)
          + log(baseline_rr_E)*S_E + log(baseline_rr_P)*S_P
log p_A = log p_C + log(rr_AC) + log(em_rr_AC)*S_E
log p_B = log p_C + log(rr_BC) + log(em_rr_BC)*S_E
```

默认乘数已经按本项目旧生成器的事件率重新校准，不再统一使用 `.67`。
z1–z3 修饰条件 RR，z4–z6 只改变基线风险。
同一个协变量组合下两项试验共享同一个 C 风险，避免公共对照不一致。
输入参数是原始总分上的风险乘数，不能直接填入旧 RR-GOP 的 alpha/beta。
三臂 16 格的风险全部检查为 `(0,1]`；无效参数直接报错，不截断概率。

独立 Bernoulli 协变量分布沿用项目原设定：BC/AD 的
`e=(.7,.3,.7,.3,.7,.3)` 固定；AC/IPD 的 q 默认按 ESS=.5 和
`direction=(-1,1,-1,1,-1,1)` 求解。可显式指定 q。
改变 ESS 不改变公共风险表。论文的 `.5/.2/.8` 人群比例没有移植进来。
每项试验仍采用随机排列的精确 1:1 分配；论文代码采用 Bernoulli(.5) 分配。

## common / rare 事件率校准（2026-09-30）

同时恢复旧 `joint_rr_gop_design(event=..., ess_ratio=.5)` 在 **AC 来源人群**和
**BC 目标人群**中的三臂平均风险，包括 MAIC 标准化的 BC 概率。
旧参数为 alpha_AC=.20、alpha_BC=.19、lambda=.30、rho=.50、beta_E=.20、
beta_P=.40；common beta0=-.26，rare beta0=-8.53。
目标值在 `data/helpers/webster_design.R` 中固定记录；新生成器不调用 OP/GOP 求根。

| 场景 | 人群 | C 风险 | A 风险 | B 风险 |
| --- | --- | ---: | ---: | ---: |
| common | AC | 40.0993% | 59.9046% | 52.9273% |
| common | BC | 41.8596% | 51.4669% | 50.2198% |
| rare | AC | 5.0131% | 7.6286% | 6.6767% |
| rare | BC | 4.8831% | 6.0809% | 5.8938% |

AC 试验实际使用 C/A，BC 试验实际使用 C/B；其余格是同一人群中的反事实风险。
这些是总体期望概率，每次生成的样本事件率有随机波动，有限样本 MAIC 也不必精确等于总体 BC 风险。

| 风险乘数/参数 | common | rare |
| --- | ---: | ---: |
| baseline_risk (b) | 0.4477141714 | 0.0459765937 |
| rr_AC | 0.8756858856 | 0.8773345625 |
| rr_BC | 1.0235503476 | 1.0244265493 |
| em_rr_AC | 1.4485437654 | 1.4637746775 |
| em_rr_BC | 1.2035936441 | 1.2104098519 |
| baseline_rr_E | 0.85 | 0.85 |
| baseline_rr_P | 1.0796973796 | 1.2278532201 |

匹配六个平均风险后仍有一个基线形状参数可选。这里固定 `baseline_rr_E=.85`，
再用独立 Bernoulli 的乘积矩求其余参数；此选择让两个场景的所有 16 格风险合法，
同时保留三个预后非 EM。`.85` 是项目设计选择，不是论文原始参数。
common 三臂格风险范围约为 C 27.50%–56.35%、A 39.21%–92.11%、B 45.83%–61.76%；
rare 为 C 2.82%–8.51%、A 4.03%–14.38%、B 4.71%–9.50%。

校准参考人群始终是默认方向、ESS=.5 的 q 和固定的 e，不随每次样本重新求参数。
因此改变 ESS、direction 或 q 后，风险表仍不变，AC 平均事件率可能改变；BC 总体风险保持固定。
显式修改风险系数也不会重新强制匹配事件率。`design$calibration` 给出目标风险、
参考风险、覆盖的参数及 `matches_reference_targets`，避免把自定义设计误称为已校准。

分数相互独立时，边际 RR 中的 `b` 和 S_P 的乘积矩抵消。
只改变非 EM 的独立分布会改变各臂风险，但不改变边际 RR。
这不意味着条件风险表或旧生成器的全部系数被复现：本次匹配的是两个参考人群的平均风险。

## 使用

```r
source("load_anchored.R", encoding="UTF-8")
source_anchored()
design <- webster_rr_design(event="common", ess_ratio=.5)
scores <- make_brm_score_spec(
  em_covariates=c("z1","z2","z3"),
  non_em_covariates=c("z4","z5","z6"))
sim <- generate_anchored_data(1000, 1000, design)
result <- run.anchored(design, 1000, 1000, seed=1,
  score_spec=scores, weight_covariates=paste0("z",1:6), run_exact=FALSE)
design$marginal_truth
result$reference_truth
```

`data/example_anchored.R` 和 `estimation/simulation_config.R` 已改用此 DGM。
旧 `joint_rr_gop_design()` 保留，显式调用它的旧脚本仍生成旧机制的数据。
旧结果不作重算或覆盖。新批量结果建议使用独立目录：

```sh
Rscript estimation/run_simulation_anchored.R estimation/simulation_config.R results_webster_rr
```

## 真值和可调场景

- 默认 A/C 和 B/C 不同；两个参考人群的 AC/BC/AB 边际 RR 均恢复旧场景水平。
- 例如 `webster_rr_design(rr_BC=.8, em_rr_BC=.9)` 给出不同 B/C 效果；
  B/C 参数不会改变 p_C 或 p_A。
- `event="rare"` 使用另一组校准参数，三臂平均风险见上表；不再只是将 common 的 b 换成 .05。
- `marginal_truth` 对 AC、BC 人群分别精确积分，提供三臂风险和 AC/BC/AB log-RR。
- `coefficient_truth` 是固定总体 h_E 尺度下的系数：
  `alpha=log(rr)+.9*log(em_rr)`、`lambda=sqrt(.63)*log(em_rr)`。
- 估计用 BC 样本均值/协方差标准化时，条件系数真值按当次分数重新计算：
  `alpha*=log(rr)+c_hat*log(em_rr)`、`lambda*=s_hat*log(em_rr)`。
  错设 EM 分数无法准确表示真实 log-RR 时不自动填写真值。
- 保留原估计目标：当前比较入口中的 marginal BC 是 BC 人群标准化的 **A/C**
  log-RR 区间。生成 B/C 数据及记录 A/B 真值，不等于已把该入口改为 A/B 间接比较。
  A/B 的 bias/coverage 不能用这些 A/C 输出直接计算。

直接 log-risk DGM 使默认 LB/LP 的风险均值形式正确，也仍能被正确 RR 加饱和
OP 的 BRM 表达。因此这是改变数据生成参数化，不意味着 BRM 被故意错设。

验证入口：`data/tests/check_webster_rr.R`、`estimation/tests/check_webster_pipeline.R`。

另有只生成数据的 [RD 机制](README_webster_RD.md) 和 `data/example_rd.R`；
它不替换默认 RR 配置，也不把当前比较流程改成 RD 分析。
