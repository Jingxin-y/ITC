# 反事实风险导出：取消概率硬筛选（2026-10-07）

范围：`cluster_all_methods.R` 导出的 BC 目标人群 A/C 边际 log RR。
对标准化人群保留协变量，分别设置 `x=0`（C）和 `x=1`（A），计算预测，
再用 MAIC 权重平均。unit 与 maic 两种拟合权重使用同一导出规则。

## 当前行为

取消 `maic_glm_target_terms()` 中按 Bernoulli FC 标签执行的 `(0,1)` 硬筛选。
原版中受影响的四个方法为 `LB-FC/RLB-FC`、`LP-BFC/RLP-BFC`；普通 `LB/RLB`、
`LP/RLP`、`LP-FC/RLP-FC` 原来没有此筛选。当前所有 GLM 方法均导出未截断的
plug-in 目标，不因逐人的反事实预测等于或大于 1 自动丢弃点估计。

仍检查：拟合成功且收敛、有效 MAIC 权重、有限且无别名的回归系数、标准化均值
有限且为正。标准化均值非有限或非正时不能计算边际 log RR，状态仍为
`standardization_failed`。方差失败时保留可独立计算的点估计，状态为 `se_failed`。
Bernoulli 事实拟合的概率约束、FC 修正方程和数值求根、标准误和 CI 规则均沿用原实现。
BRM 的 RR/OP 概率映射、CMH 的加权分层计数和条件系数导出不变。

没有截断预测、删除标准化个体或改变权重。没有新增逐人越界率的自动诊断字段；
这次修改仅取消硬筛选。可计算的点估计不表示每个预测是合法概率，边际风险均值
也可能大于 1。概率有效性可另做诊断或敏感性分析，不混同于拟合是否收敛。

例如只观察到 `(x,z)=(0,0),(0,1),(1,0)`，风险分别为 `.2,.8,.6`。
无交互 log-link 模型可以收敛并精确拟合，但未观察到的 `(1,1)` 预测为 `2.4`。
BC 中 `z=1` 占 1% 时，标准化均值为 `.206,.618`，边际 log RR 为 `log(3)`。
当前规则保留这个可计算的估计；反事实预测越界仍是概率有效性问题。

## 检查点和集群更新

`settings.version=10`，保存
`counterfactual_risk_check="finite_positive_standardized_means"`。
版本 8 的 FC 硬筛选和版本 9 的统一 Bernoulli 硬筛选不能作为版本 10 续跑。
新运行使用新结果目录，历史图片、CSV、检查点和代码快照不会自动修改。
仅为本次运行逻辑修改，集群需要更新 `estimation/cluster_all_methods.R`；
其余改动是核验脚本和说明，不是运行依赖。

核验：`Rscript estimation/tests/check_counterfactual_risks.R` 检查真实拟合的反事实
预测越界后仍保留边际点/SE/CI，unit/maic 导出，概率等于 1 的边界，合法预测结果不变，
数值不可计算仍失败、方差失败恢复点估计、Bernoulli 事实域约束及 FC 未收敛仍失败。
