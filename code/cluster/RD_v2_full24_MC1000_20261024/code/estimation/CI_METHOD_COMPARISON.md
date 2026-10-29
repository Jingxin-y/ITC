# 逐系数 BC 与线性目标 BC：源码比较

核对日期：2026-09-19。逐系数版本为 [brmplus/compare/CI_exact_adapt_para.R](https://github.com/hta-pharma/brmplus/blob/3b628a6d32dc7fe0f31fbd3ab181aecd3866f982/compare/CI_exact_adapt_para.R)，commit `3b628a6d32dc7fe0f31fbd3ab181aecd3866f982`。本地保存的该文件与上游内容一致。线性目标版本为本目录下 `helpers/CI_linear_target.R`，本次保留其数值实现不变。

本文分析的是代码实际执行的统计构造及数值实现，不是有限样本覆盖率证明。

## 1. 本次文件整理和调用关系

- `helpers/CI_exact_adapt_para.R`：上游原版，定义 `exact()` 和并行 worker。
- `helpers/CI_linear_target.R`：保留，定义 `exact_linear_target()`。
- 删除旧串行移植、旧 RLP 文件名占位和 4 小时限时包装；时限直接传 `time.limit.sec`。
- `source_anchored()` 默认同时加载两种函数。当前主流程默认显式调用 `exact_marginal_target()`。本文逐系数调用链仅适用于 `bc_target="coefficient"`；线性目标仍需显式调用。
- 上游 `exact()` 不接受旧移植版 `B.*` 参数。它的 bootstrap 次数固定为 100/200/200/500；线性版通过 `nsim.*` 参数调整次数。

## 2. 推断目标：单个系数是线性目标的特例

两者使用同一类结局模型：

\[
g\{E(Y\mid X=1,Z),E(Y\mid X=0,Z)\}=v_a(Z)^T\alpha,
\qquad \log OP(Z)=v_b(Z)^T\beta.
\]

其中 RR 模型的效应链接是 `log(RR)`，RD 模型是 `atanh(RD)`；beta 是 odds-product 的 nuisance 参数。

| 项目 | 上游 `exact()` | `exact_linear_target()` |
|---|---|---|
| 目标 | 每一个 `alpha[j]` | 一个用户指定的 `psi = c^T alpha` |
| 一次调用 | 循环计算 `pa` 个系数的区间和 P 值 | 计算一个 contrast 的区间和 P 值 |
| 原假设 | `alpha[j] = 0` | `psi = psi.null`，默认 0 |
| 置信水平 | 固定 95%，代码阈值为 0.05 | 可配置 `level`，默认 0.95 |
| 其他 alpha 系数 | 给每个目标作 profile 时重新拟合 | 给定 `c^T alpha` 后在剩余方向重新拟合 |
| beta | 重新拟合 | 重新拟合 |
| 输出尺度 | alpha 的链接尺度 | 链接尺度，并额外给出 `exp(psi)` 或 `tanh(psi)` |

“逐系数”不表示模型只有一个 alpha，也不表示忽略 alpha 之间的相关性：两者都在完整似然中联合处理 nuisance 参数。逐系数版本输出的是各个系数的点态区间，不是联合置信区域。

若 `c=e_j`（第 j 个位置为 1，其余为 0），则约束 `c^T alpha=psi` 就是 `alpha[j]=psi`。在相同可行域、充分准确优化、相同检验阈值和相同模拟分布下，两种构造应对应同一个目标；目前实现的可行域、停止判据、调用顺序和浮点阈值存在差异，因此不能要求当前代码逐次输出相同。

例如 `log RR(z)=alpha_0+z*alpha_1`：

- `c=(1,0)`：`z=0` 时的条件 log RR。
- `c=(1,0.3)`：`z=0.3` 时的条件 log RR。
- `c=(0,1)`：两个相差 1 单位的 z 取值之间的 log RR 差；取 exp 后是两个条件 RR 的比值。

不能把两个单系数区间的下限相加、上限相加来代替 `alpha_0+0.3*alpha_1` 的 BC 区间；应直接约束并反演该组合。

## 3. 共同的统计构造

记加权负对数似然为

\[
Q(\alpha,\beta;y)
=-\sum_i w_i\{y_i\log p_i+(1-y_i)\log(1-p_i)\}.
\]

两者的 `p_i` 都来自 `getProbRR()` 或 `getProbRD()`，并裁剪到 `[1e-12,1-1e-12]`。给定目标函数 `h(alpha)`，定义其 profile 值：

\[
Q_t(y)=\min_{h(\alpha)=t,\;\beta} Q(\alpha,\beta;y).
\]

逐系数版取 `h(alpha)=alpha[j]`，线性版取 `h(alpha)=c^T alpha`。实际代码通过有界交替优化近似求这个最小值。

令 \(\hat t=h(\alpha_{input})\)，其中 `alpha_input` 来自传入的 `pars`；它可能是 MLE、Firth 估计或外部修正估计，虽然变量名含 `.ml`，不能因此把所有输入都称为无约束 MLE。

对每个候选 `t`，两者都执行：

1. 在观察数据中拟合 `h(alpha)=t`，得到候选概率。
2. 计算 \(T_{obs}(t)=2\{Q_t(y)-Q_{\hat t}(y)\}\)。
3. 从候选受约束拟合生成独立 Bernoulli 结局 `y*`，保持 X、Z 和权重不变。
4. 对每个 `y*` 分别拟合目标为 \(t\) 和 \(\hat t\) 的约束，计算 \(T^*(t)=2\{Q_t(y^*)-Q_{\hat t}(y^*)\}\)。
5. 用以下经验 acceptability 判断候选值是否接受。

关键点：每个 bootstrap 样本中，参照目标值 \(\hat t\) 仍固定为原始输入点估计对应的值，没有在模拟样本上重新寻找一个无约束目标估计。两个受约束拟合的差也没有被强制为非负。因此，它不同于通常“受约束拟合减去该样本无约束 MLE 拟合”的标准 profile LRT。

所有新的候选拟合和 bootstrap 拟合都从原始输入参数开始；缓存减少同一候选的重复观察数据拟合，但不把相邻候选的解当作 warm start。

### Acceptability 的具体公式

对模拟值 `T_1*,...,T_B*`，定义

\[
F_B(t)=B^{-1}\sum_b I(T_b^*\le t),\quad
S_B(t)=B^{-1}\sum_b I(T_b^*\ge t),\quad
m_B(t)=\min\{F_B(t),S_B(t)\}.
\]

两份代码均使用

\[
A_B(t)=B^{-1}\sum_b I\{m_B(T_b^*)\le m_B(T_{obs})\}.
\]

这是先取左右尾概率的较小值、再按其模拟分布校准的 Blaker-type acceptability。区间由接受集合的端点搜索产生。P 值是将目标设为原假设值后再次计算同一个 acceptability。

有限 B 下 acceptability 按 `1/B` 的步长变化，并有 Monte Carlo 随机性；函数名 `exact` 本身不构成有限样本精确覆盖率的证明，尤其 nuisance 参数由数据估计且拟合、搜索均有数值容差。

## 4. 最重要的数值差异：约束和优化可行域

### 逐系数版

每次直接将 `alpha[j]` 设为候选值，在目标函数内部再次固定该分量；`optim()` 仍收到长度 `pa` 的 alpha 向量，因此其中一个方向实际是平的。交替进行：

- alpha 优化：L-BFGS-B，逐坐标界限 `[-8,8]`。
- beta 优化：L-BFGS-B，逐坐标界限 `[-10,10]`。

目标候选本身也由区间搜索限制在 `[-8,8]`。

### 线性版

用 QR 分解得到 `c^T` 的零空间正交基 N，并写成

\[
\alpha=\frac{c}{c^Tc}\psi+N\gamma,\qquad c^TN=0.
\]

所以任意 gamma 都满足 `c^T alpha=psi`，只需优化 `pa-1` 个 gamma 和 beta。`pa=1` 时不再进行 alpha/gamma 优化。

初始化时，线性版把原输入 alpha 正交投影到候选约束平面；逐系数版只替换待检验的一个坐标。单位 contrast 时这两种初始化描述相符，一般 contrast 则可能同时改变多个 alpha 坐标。

但当前默认界限是

\[
|\gamma_k|\le 8\sqrt{pa},\qquad |\beta_k|\le10.
\]

**这个 gamma 盒约束不等价于原来的逐坐标 `|alpha_k|<=8`。** 代码也没有在重建 alpha 后强制检查该界限。`alpha.bound` 主要用于确定默认 gamma 界限和目标搜索界限。

简单例子：`pa=2,c=(1,0)` 时，自由方向是第二个 alpha（可能相差一个符号）；线性版允许它到约 `±11.31`，上游只允许 `±8`。在边界附近，两者可能得到不同 profile 最优值。这种差异即使增大 bootstrap 次数也不会自动消失。

若要严格统一可行域，应在 gamma 空间施加 `-8 <= c*psi/(c^T c)+N*gamma <= 8` 的线性不等式，或明确决定两种实现都采用同一个新可行域。单纯设置 `gamma.bound=8` 不能解决一般 contrast 的等价性问题。

### 收敛判据也不同

上游使用整个 alpha 向量的相对平方变化：

`sum((alpha_new-alpha_old)^2) / sum(alpha_new^2 + thres)`。

固定目标分量也进入分母。线性版改为 gamma 的相对平方变化，并取其与 beta 变化的最大值。二者即使描述相同约束，也可能在不同迭代停止。L-BFGS-B 的有限差分及坐标变化也会产生数值差异。

两者都记录优化失败，并在错误/非有限参数等情况下回退到起始参数，没有自动重试或保证全局最优。profile 的 `converged` 主要反映交替参数变化达到阈值；它不能代替检查 `optimizer.failures`。

## 5. 区间搜索和标准误的作用

| 项目 | 逐系数版 | 线性版 |
|---|---|---|
| 中心 | 输入 `alpha[j]` | 输入 `c^T alpha` |
| 搜索标准误 | `se[j]` | 优先用户给的 `psi.se`；否则 `sqrt(c^T vcov.alpha c)`；最后才用忽略协方差的 SE 后备 |
| 默认半宽 | n<300 用 8SE，否则 4SE | 相同默认，还允许 `search.mult` |
| 目标界限 | `[-8,8]` | `[-8*sum(abs(c)),8*sum(abs(c))]`，可配置 `alpha.bound` |
| 默认网格比例 | `0,.04,.08,.16,.32,.64,1` | 相同 |
| 粗搜索/确认/细化/P 值 B | 100/200/200/500，内部常量 | 默认相同，`nsim.*` 可配置 |
| 搜索确认 | 粗拒绝点复查；一般需再一个向外拒绝点，最外端允许单个确认拒绝点 | 相同规则 |
| 二分搜索 | 默认宽度 <=0.001 或最多15步 | 默认相同 |
| 端点取值 | 最后的接受侧 inner 点 | 相同 |
| 搜索顺序 | 先所有系数上限，再所有下限，再所有 P 值 | 单个目标先下限、再上限、再 P 值 |

SE/协方差用来确定**有限搜索范围**，区间本身由 acceptability 反演。故不是把 SE 直接代入 Wald 公式；但若 SE 太小导致拒绝边界落在搜索范围外，仍会影响最终输出。代码此时返回搜索边界并标记不收敛，而不会无限扩大范围。

对于线性目标 `alpha_0+z0*alpha_1`，正确的搜索方差含 `2*z0*Cov(alpha_0,alpha_1)`。只提供边际 SE 的后备会省略这项，但 profile 拟合本身仍保留参数间的似然依赖关系。

两者只搜索围绕输入估计的两侧端点；这个数值过程不是对可能不连续或非单调的整个接受集合的穷尽刻画。

## 6. 已复现的默认阈值差异

上游直接比较 `A >= 0.05`。线性版先计算 `alpha.level <- 1-level`，默认比较 `A >= 1-0.95`。

R 双精度下：

```text
0.05      = 0.050000000000000003
1 - 0.95  = 0.050000000000000044
0.05 >= 0.05      -> TRUE
0.05 >= 1 - 0.95  -> FALSE
```

当 B=100 或200时，经验 acceptability 可以恰好为0.05；这时上游接受、当前线性版拒绝。该差异已经用实际 R 比较及 5/100、10/200 的 `mean()` 复现，不能仅归因于随机种子不同。

此外，线性文件开头第18行写“acceptability >= level”，但实际执行的是 `>= 1-level`。该注释应改为 `>= 1-level`，而数值比较需要明确浮点容差处理。这次只记录问题，没有修改线性方法的数值逻辑或其原始注释。

## 7. 随机数、并行、输入检查和返回值

- 上游可串行、Windows PSOCK 或 Linux/macOS fork；默认串行。线性版只有串行。
- 上游先生成一批 bootstrap 结局，再拟合；线性版生成一个就拟合一个。在拟合及概率助手不使用随机数的前提下，这一点本身不会改变对应批次的抽样序列。
- 但上下限计算顺序不同，且上游会计算所有系数；因此即使 `set.seed()` 一致，同一候选通常拿不到同一批随机数。不能以同种子端点完全相等作为现有两种实现的等价要求。
- 线性版增加 `seed` 参数，调用 `set.seed()` 后会改变调用环境的全局随机数状态；上游由调用者设置种子。
- 线性版检查 contrast 非零、矩阵维度、param 属于 RR/RD、level 范围等。上游对时间上限、二分步数、`grid.mult>1`、二分容差有显式检查；线性版缺少对应的完整检查。例如 `grid.mult=1` 会令构造网格的 while 循环不能前进，且该循环内部不检查时限。
- 上游超时使用 `exact_timeout` 条件；线性版使用普通错误。仅按条件类捕获超时的外部代码不能直接通用。
- 上游输出长度为pa的 `low/up/p` 和汇总 diagnostics；线性版输出单个 `estimate/low/up/p`、contrast、变换后的 effect 及约束误差等诊断。接口不能直接替换。

## 8. 条件效应、边际效应和 RD 变换的解释

`c=v_a(z0)` 时，`exp(c^T alpha)` 或 `tanh(c^T alpha)` 正好是指定协变量取值 z0 下的条件 RR/RD。

若把 c 设为目标人群平均设计向量：

\[
\exp\{E[v_a(Z)]^T\alpha\}
\]

对应平均 log RR 的指数（在该线性模型下），一般不等于边际 RR：

\[
RR_{marginal}=\frac{E[p_1(Z)]}{E[p_0(Z)]}.
\]

RD 同理，`tanh(E[v_a]^T alpha)` 一般不等于 `E[p1(Z)-p0(Z)]`。目标人群边际效应通常是 alpha、beta 和目标协变量分布的非线性函数，当前线性约束函数不能自动完成这种推断。

对任意 contrast，自动输出的 `effect.*` 只能先理解为链接目标的单调变换。尤其 RD：若 `psi=atanh(RD(z1))-atanh(RD(z0))`，则 `tanh(psi)` **不是** `RD(z1)-RD(z0)`。必须先定义 estimand，再解释这些字段。

## 9. 两者共有的推断范围

- bootstrap 固定治疗分组、协变量、MAIC 权重；只模拟结局，没有重新估计权重或抽样 AD 汇总量。
- 线性版 contrast 也被当作固定值。若它来自估计的目标均值，当前方法不传播 contrast 本身的估计误差。
- 输入 joint sandwich SE/协方差只改变搜索宽度，不会自动把权重或 AD 均值的不确定性加入 bootstrap。
- 两者 profile 的目标函数都没有 Firth 惩罚。传入 Firth 参数会改变输入锚点、初值和搜索范围，但不会把 bootstrap 内部拟合转换为 Firth 拟合。
- 搜索端点收敛、优化器成功和统计覆盖率是不同问题，不能用运行完成替代这些检查。

## 10. 使用方式和后续验证

```r
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

# 逐系数：保留上游默认模拟次数、算法和串行设置。
ci.coefficients <- exact(
  param = "RR", y = y, x = x, va = va, vb = vb, weight = weight,
  max.step = 50, pars = pars, se = se, pa = ncol(va), pb = ncol(vb)
)

# 条件效应：假设 va 的两列按顺序为截距和 z。
ci.at.z <- exact_linear_target(
  param = "RR", y = y, x = x, va = va, vb = vb, weight = weight,
  contrast = c(1, 0.3), max.step = 50, pars = pars,
  pa = ncol(va), pb = ncol(vb), vcov.alpha = V.alpha, seed = 1
)
```

逐系数方法现直接使用上游文件。线性版保留了相同的校准构造，但正式用于方法比较前，优先统一/修正：阈值的浮点边界、alpha 可行域、优化停止规则、输入检查，以及所需 effect 字段的 estimand 解释。验证可先用单位 contrast 比较相同候选处的 profile 值和 acceptability，再检查端点；只比较一次随机区间，无法区分优化差异和 Monte Carlo 波动。
