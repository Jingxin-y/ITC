# 独立 RD 生成器：三个 EM、三个非 EM

`webster_rd_design()` 使用直接加法风险模型，借鉴
[Webster-Clark et al. (2024)](https://doi.org/10.1093/aje/kwae048) 的 RD 思路。
这是保留本项目三 EM / 三非 EM、两试验、分数方向和 ESS 设定的扩展，
不是论文原始参数的复现。默认示例和批量估计仍用 RR 生成器；本次没有改动 RR 估计实现。

```text
S_E = (1-z1)+z2+(1-z3)
S_P = z4+(1-z5)+z6
p_C = b + beta_E*S_E + beta_P*S_P
p_A = p_C + rd_AC + em_rd_AC*S_E
p_B = p_C + rd_BC + em_rd_BC*S_E
```

所以 `RD_AC=rd_AC+em_rd_AC*S_E`、`RD_BC=rd_BC+em_rd_BC*S_E`。
z1–z3 是 RD 的 EM，z4–z6 改变基线风险但在风险差中抵消。
边际 RD 为 `rd+em_rd*E[S_E]`：只改变 S_P 分布并保持 S_E 分布，边际 RD 也不变，
不要求 S_E/S_P 独立。RD 的非 EM 不自动成为 RR 的非 EM，这些数据应按其目标尺度解释。

默认 common/rare 同时匹配原 RR-GOP 在 AC、BC 人群中的三臂平均事件率，
参见 [RR 机制中的事件率表](README_webster_RR.md)。校准参考 ESS=.5/default direction，
六变量仍独立生成。固定基线形状 `beta_P=.02` (common) / `.002` (rare)，
令 mu_AC、mu_BC 为 S_E/S_P 的两个人群均值，然后求：

```text
beta_E = (target_C_AC-target_C_BC)/(mu_AC-mu_BC) - beta_P
b = target_C_BC - (beta_E+beta_P)*mu_BC
em_rd = (target_RD_ACpopulation-target_RD_BCpopulation)/(mu_AC-mu_BC)
rd = target_RD_BCpopulation - em_rd*mu_BC
```

其中 `target_RD` 分别取该人群 A−C 或 B−C。这些选择保留三个预后非 EM，
且 common 三臂格风险范围为 C 28.81%–50.88%、A 36.98%–85.27%、B 45.57%–61.07%；
rare 为 C 4.66%–5.40%、A 3.42%–12.28%、B 4.55%–9.03%。
所有 16 格三臂概率必须严格在 `(0,1)` 内，无效参数报错，不截断或取 logistic，
避免破坏加法 RD 关系。单个 RD 可为负；默认三个变量的 RD 修饰系数均非零。

```r
source("load_anchored.R", encoding="UTF-8")
source_anchored_data()
design <- webster_rd_design(event="rare", ess_ratio=.5)
set.seed(1)
sim <- generate_anchored_data(1000,1000,design)
design$marginal_truth[,c("p_C","p_A","p_B","RD_AC","RD_BC","RD_AB")]
design$rd_coefficient_truth
design$calibration
```

共享抽样器照常返回 `data.IPD`、`data.AD.or` 和 `data.AD`，两项试验各为精确 1:1 随机分配。
`rd_coefficient_truth` 只表示原始 RD 在固定 h_E 单位上的截距/斜率，不是现有 RR 模型的 alpha。
风险表与 `marginal_truth` 中保留实际 log-RR 作描述，没有将其声明为线性 log-RR 真值。

改变 ESS/q 时不重新校准风险参数；显式覆盖系数时，其余参数保留默认值，
`calibration$matches_reference_targets` 会反映参考事件率是否仍匹配。
估计用权重、MAIC 标准化和 RR 比较方法均按现有设置运行，没有新增 RD 比较方法。

完整示例：`data/example_rd.R`；检查：`Rscript data/tests/check_webster_rd.R`。
