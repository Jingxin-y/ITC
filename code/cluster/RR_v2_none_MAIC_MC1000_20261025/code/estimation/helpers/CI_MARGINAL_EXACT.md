# Marginal logRR 的 bootstrap 检验反演区间

入口：`CI_marginal_target.R` 中的 `exact_marginal_target()`。
`source_anchored()` 会加载它；也可以直接 source 此文件，求解器仅依赖 base R。

目标是

```
psi = log(sum(v * p_A) / sum(v * p_C))
```

其中 `p_A,p_C` 来自两臂 RR + OP 模型。函数同时优化 RR 与 OP 参数，
不把 marginal logRR 当作截距或 alpha 的线性对比。

## 调用已有 BRM / FC 拟合

```r
source("load_anchored.R", encoding="UTF-8")
source_anchored()

# fit 是已有全 1 权重的 brm 或 brm-FC 拟合；dat 与设计矩阵的行顺序一致。
ci <- exact_marginal_target(
  y = dat$y, x = dat$x,
  va = fit$va, vb = fit$vb,
  pars = fit$point.est,
  weight = rep(1, length(dat$y)),
  standardization.weight = maic_weights,
  vcov.pars = fit$cov,
  seed = 100001L
)
ci[c("estimate", "low", "up", "p", "status")]
ci[c("effect.estimate", "effect.low", "effect.up")]
ci$diagnostics$counts
```

`weight` 是拟合权重，`standardization.weight` 是目标标准化权重。
当前 BRM/FC 默认全 1 权重拟合，区间应传 `weight=rep(1,n)`，
同时显式传入 MAIC 标准化权重。这两个参数不可混淆。
`weight` 用于观察数据及 bootstrap 数据的 nuisance likelihood；
`standardization.weight` 用于 marginal 风险平均和每次截距约束求根。
底层函数省略标准化权重时默认使用拟合权重，因此 unit 模式下不可省略
MAIC 目标权重。若用完整主流程结果，分别传
`result$model_weights$weight` 和 `result$weights$weight`。

若使用单独的标准化人群，可同时提供 `va.target`、`vb.target` 和
`standardization.weight`；设计列必须与原拟合含义、顺序一致。
仅有 BC 汇总量时，不能把协变量均值作为“代表个体”代替联合分布。

可运行的真实项目示例：

```r
source("estimation/example_marginal_exact.R", encoding="UTF-8")
ans <- example_marginal_exact(method="brm")
# 或 method="brm-FC"
ans$ci[c("estimate","low","up","p","status")]
```

示例默认在估计阶段以 z1-z3 求 MAIC 目标权重，BRM/FC 用全 1 权重拟合，
正确使用 h_E/h_P、保留原数据生成；区间使用与拟合相同的全 1 likelihood 权重。
可传入 `weight_covariates=paste0("z",1:6)` 改用全变量权重。
主比较及 simulate_anchored() 已接入本函数：`run_exact=TRUE` 默认使用
`bc_target="marginal"`。`result$marginal_estimates` 显示 marginal logRR 结果；
`result$bc` 保存原始返回值和诊断。旧逐系数区间需显式指定
`bc_target="coefficient"`。旧批量 grid 的历史脚本及已有结果没有改写。

主流程默认 `model_weighting="unit"`：拟合权重为全 1，标准化权重为本次
MAIC 权重；BRM/FC 复用各自拟合的 `va`、`vb` 和 `point.est`。
独立 `get.brm.bc()` 也默认 unit，需提供 unit 拟合的输入模型。
显式 `model_weighting="maic"` 保留历史加权拟合方式。
Bayes 修正仍仅用于无 EM 情形。
BC 行标记 `component="marginal"`、`parameter="log_RR"`；原 Wald 系数行保留。
bootstrap 不返回标准误，因此 BC 行的 `se=NA`、`se_status="not_provided"`，
不会把 alpha SE 或仅用于搜索的 delta SE 当作 marginal bootstrap SE。
`bc_control` 接收本函数的搜索/模拟参数，不允许覆盖数据、模型或标准化权重。
默认 `optim.maxit=200`；本函数没有 `parallel.bootstrap` 接口。

批量 CSV 另输出 `anchored_joint_marginal_exact.csv`。汇总包含 marginal 行，
但不把 alpha 系数真值用于 marginal 覆盖率；需要显式提供
`method/component="marginal"/term="marginal log RR"/truth` 对应的目标真值。

## 计算方法

每个候选 `psi0` 下，以其他 RR 系数和全部 OP 系数为 nuisance 参数。
通过一维求根确定 RR 截距，令 marginal logRR 严格等于 `psi0`。
OP 参数每次变化时也重新求截距，不能只在优化开始时施加一次约束。
默认 unit 流程在 MAIC 目标约束下优化全 1 权重 likelihood；
bootstrap 两次拟合沿用相同规则，并不把候选 psi 的 likelihood 改为 MAIC 加权。

给定非截距 RR 线性预测 `offset`，截距根的解析括区间为
`[psi0-max(offset), psi0-min(offset)]`，只使用标准化权重为正的个体。
拟合使用隐函数求导得到的解析约束梯度。RR 截距必须对应全为 1 的设计列，
默认 `intercept=1`。2026-09-27 按用户要求恢复为 2026-09-22 结果对应的
原版 Exact 源文件（MD5 `18e26a7b958a84ae0421fec6c8cc7341`）。
2026-10-02 按用户要求，仅调整 nuisance 起点和默认边界：使用有上下界的
`L-BFGS-B` 联合优化，非截距 alpha 和全部 beta 的默认范围均为 `[-15,15]`。
越界的优化起点移到对应边界内侧；原始 `pars` 保留，仍用于计算观察数据的
marginal 点估计。RR 截距继续通过 marginal 约束求解，不受上述边界限制。
显式指定的 `alpha.bound`、`beta.bound` 仍生效；不自动扩展边界，
也不接受 `expand.bounds`。OP 编码和两套权重的用途保持原有规则。
默认 `optim.reltol=1e-10`，通过 `factr=optim.reltol/.Machine$double.eps`
传递；同时设置 `pgtol=score.tol/10`。第一次优化的上限仍为 200。
收敛要求优化器返回 code 0、投影 nuisance score 不超过 `score.tol=1e-4`，
以及 marginal 约束残差不超过 `constraint.tol=1e-7`。

2026-10-04 增加失败拟合重试：`optim.retries=2`，即最多三次求解。
只对未通过检查的拟合重试；保留同一份 bootstrap 结局、同一候选 psi、
同一组拟合及标准化权重。从上次 nuisance 参数（发生错误时从最后有效
目标值评价位置）继续 `L-BFGS-B`。每次重试上限默认
`optim.retry.maxit=max(1000,optim.maxit)`，目标值停止容差依次为
`optim.reltol/100`、`optim.reltol/10000`，通常下限为机器精度；
若显式指定了更小的原容差，则保留该容差作为下限，不放宽它。
优化器内部 `pgtol` 同样依次收紧 100 倍、10000 倍，防止极接近边界时
其投影步长很小而提前停止。
投影 score、约束残差的通过阈值及参数边界保持原值。
重试不抽新数据、不丢弃失败重复，也不增加 bootstrap 样本数。
`bc_control=list(optim.retries=0L)` 可以关闭重试，用于旧行为对照。
重试耗时计入原来的 `time.limit.sec`；超时继续向上返回。

`settings` 记录重试设置、边界及阈值；拟合保留
`max.projected.score`、`nuisance.boundary`、`optimizer.attempts`，
重试过的拟合另有 `attempts` 和 `retry.recovered`。
辅助参数达到边界与 marginal 区间端点达到 `psi.range` 是两种不同状态。
区间的 `boundary_truncated` 标记继续保留。

在候选约束拟合下生成 Bernoulli 结局。每份 bootstrap 数据分别拟合候选
`psi0` 和固定的观察数据中心 `psi.hat`，用两者的加权负 log likelihood 差
计算统计量，再按原 `exact()` 的 Blaker-type acceptability 规则反演。
这不是每份 bootstrap 都重新寻找无约束 MLE 的常规 LR bootstrap。

FC 参数可以作为输入中心，沿用现有 exact 的做法；nuisance 拟合使用普通
加权 likelihood，不声称实现了 FC 惩罚 likelihood。

概率变换采用 RR/OP 方程的稳定二次根，与现有 `getProbRR()` 在正常范围
一致，但不采用其超过 +/-12 时的边界近似。原函数没有被修改。

## 返回值与状态

- `estimate,low,up`：marginal **logRR**；`effect.*`：对应 RR。
- `p`：对 `psi.null`（默认 0）的 bootstrap acceptability。
- `interval.available`：是否获得两个可保留端点，包括显式边界截断。
- `status="boundary_truncated"`：在用户指定/默认搜索界限仍接受，返回该界限，
  并记录截断；未确认区间真正到此结束，不自动扩展搜索范围。
- `interval_failed`：至少一端出现拟合或搜索失败。失败端点为 NA，
  不使用末次未收敛参数伪装成有效区间。
- `diagnostics$failed.fits`：最多保存前 25 次重试后仍失败的拟合、原因及
  结局向量 `response`；`diagnostics$retried.fits` 同样最多保存前 25 次
  重试过的拟合及结局。`counts$retried.profile.fits` 记录重试拟合数，
  `counts$recovered.profile.fits` 记录恢复数，`counts$optimizer.attempts` 记录
  实际优化次数；`counts` 原来的失败次数只统计重试后仍失败的拟合。
  bootstrap 中任何拟合最终失败都会使该候选的
  acceptability 无效，不能删除失败重复后重新计算比例。
- `diagnostics$candidates`：观察数据约束拟合及其约束残差。
- `diagnostics$evaluations`：各候选 acceptability、模拟次数及失败状态。
- `diagnostics$lower/upper`：端点括区间、最后接受点及搜索状态。

`psi.grid` 可以额外计算指定候选的 acceptability；
`compute.interval=FALSE`、`compute.p.value=FALSE` 可跳过相应计算。
当前端点算法沿用由中心向外搜索及二分，不保证识别非连通接受集合的全部分支。

## 精度、迭代和统计口径

默认 `optim.maxit=200`；兼容参数 `max.step` 仅作为 `optim.maxit` 的默认值。
新实现联合优化全部 nuisance 参数，没有旧代码的外层交替循环。
`root.maxit=200` 是截距一维求根上限，与 optim 迭代次数不同。
普通 BRM 点估计是另一个求解器，其主流程默认
`fit_control$brm_max_step=1000`（原为 200），内层 `brm_optim_maxit=200`、
`brm_thres=1e-6` 和 `brm_optim_reltol=1e-6` 保持原值。
FC 与 GLM 的步数上限保持原值；增加上限只允许更多迭代，仍按原条件判断收敛。

默认 bootstrap 次数沿用原线性目标版本：coarse=100、confirm=200、
refine=200、p.value=500。它们是起始配置，不表示端点达到高精度；在 5%
尾部做正式比较时应增加模拟次数，并检查不同种子和更宽搜索范围的稳定性。
`thres.dicho` 仅控制数值括区间宽度，不能消除 Monte Carlo 误差。

`psi.se` 或完整的 `vcov.pars` 仅确定初始搜索宽度，不用于形成 Wald 端点。
也可直接指定包含点估计的 `psi.range=c(lower,upper)`，此时不必提供 SE。
不能仅传入 alpha 的协方差而遗漏 OP 参数的不确定性。

本版本固定协变量、治疗分配、拟合权重和标准化权重，只生成结局。
它与包含权重/BC 汇总量不确定性的 joint Wald 不是同一不确定性口径。
函数名中的 exact 延续现有命名；实际方法是 plug-in nuisance 参数下的
参数 bootstrap 检验反演，不承诺严格有限样本精确覆盖。

## 核验

```r
# 在 anchored_v1 根目录执行
# Rscript estimation/tests/check_marginal_exact.R
# Rscript estimation/tests/check_ci_weight_separation.R
# Rscript estimation/tests/check_marginal_retry.R
```

核验包括风险方程、全参数梯度、隐函数约束梯度、独立固定系数 profile 对照、
完整区间反演、权重缩放和随机种子、失败保存、截断保留与超时传播。
测试的小 bootstrap 次数用于验证程序行为，不用于验证 95% 覆盖率。
