# Rare 数值处理更新（2026-09-23）

2026-09-27 恢复更新：当前 marginal Exact 已恢复到 2026-09-22 原版，
使用有界 L-BFGS-B、alpha.bound=8、beta.bound=10、optim.reltol=1e-10、
score.tol=1e-4、optim.maxit=200。下文的自动扩展边界方案及此前的
无界 BFGS 均为历史记录；当前不接受 expand.bounds。普通 BRM 内层显式使用
`reltol=1e-6`，外层阈值仍为 `1e-6`；joint Wald 的 Ag/H 已使用
`MASS::ginv()`。当前规则见 `FAILURE_HANDLING_AND_SCORES.md` 和
`helpers/CI_MARGINAL_EXACT.md`。

本次只更新估计代码；现有两批 rare 的原始文件、汇总和图均保持不变。
已有结果不会因为修改代码自动获得新区间或新的有效次数。

## 点估计与 SE 分开保留

旧的 cluster marginal 包装函数在 `maic_delta()` 中先计算点估计，随后
计算联合 sandwich 方差。`Singular joint bread.` 的异常使整个返回值丢失，
外层预置的点估计 NA 没有被覆盖。

现在 GLM marginal 标准化的点估计可以独立计算。若点估计存在、方差失败，
保留 estimate/risk_C/risk_A，标记 `se_failed`，记录原错误；SE 和区间仍为 NA。
已核验两种权重下各一个 RLP 和 RLB 的真实失败案例，恢复的点估计均与
相同 seed 的对应 LP/LB 点估计一致。现有结果未被追溯填补。

有限的大 SE 可以作为超出图示范围的点计数，仍参与均值。NA/NaN 没有方向，
应单列缺失数量。计算溢出的 Inf 应单列非有限数量，不能据此制造无限宽区间；
奇异矩阵也不等价于标准误等于正无穷。统计模型真实返回的单侧/无界区间属于
另一种情况，需要保留其明确状态与端点定义。

## Exact 数值范围

核对两批全部 9,000 个数据重复的系数 CSV，以及 24,790 个有效 BRM 类拟合：

- 非截距 alpha 绝对值最大 5.3638，保留基础 `alpha.bound=8`。
- 各组 beta 绝对值最大值的第 99 百分位最高 16.3003；基础 `beta.bound` 从 10 改成 20。
- beta 极端最大 139.2939；103 个拟合仍超过 20，不把所有拟合统一设成极宽范围。
- `expand.bounds=TRUE` 时，对输入系数绝对值达到基础边界 90% 的坐标，
  使用 `max(基础边界, 1.25*abs(输入系数)+2)`，为起点留出余量。
- 有效边界从观察数据输入决定一次，之后所有候选和 bootstrap 重复保持相同。
  请求边界、有效边界、扩展标记及起点均保存到 diagnostics。
- 通过 `expand.bounds=FALSE, alpha.bound=8, beta.bound=10` 可使用旧规则。

新规则对 3 个旧超界拟合完成了小预算 bootstrap 区间试跑，包含 beta=139.2939
的极端案例，均返回有效区间。此项验证只支持程序运行和状态处理正确，不能证明
95% coverage 改善。扩大范围仍可能遇到拟合失败或不稳定。

有效的 `boundary_truncated` 区间继续计入 coverage。它表示已计算到搜索界限，
在该处仍被接受；它不证明真正区间止于该界限。起始辅助参数超界、bootstrap
失败或端点为 NA 属于不同情况，不能当成已经得到一个宽区间。

## 收敛阈值

在每个有失败的情形/方法分层，固定选前两个失败 seed，共检查 71 个未收敛拟合：
50 个普通 BRM 的内层优化仍未成功，21 个 FC 的末次 scaled score 仍大于 1e-3。
这不是随机抽样，也不是重新拟合实验；它说明这些保存的末次结果并不能只靠
把阈值改为 1e-3 就变成收敛结果。重新运行后的停点可能改变，需做敏感性比较。

新增可配置的 `brm_thres` 与 `fc_thres`，默认保持 1e-6。示例：

```r
fit_control = list(brm_thres=1e-3, fc_thres=1e-3)
```

可传给 `run.anchored()`、`simulate_anchored()` 或 `maic_run_seed()`。
BRM 的阈值衡量平方相对参数变化，FC 衡量标准化调整得分，两者含义不同。
这不会放宽内层 optim 收敛、Exact score、约束或求根阈值。内层 optim 上限仍为 200。

## 验证与诊断材料

新增 `estimation/tests/check_rare_boundaries.R`，验证扩展范围、记录与回退、
阈值传递，以及奇异方差不能变成伪区间。现有 marginal exact、marginal pipeline、
迭代保存测试也用于回归验证。

参数分布、71 个末次拟合诊断、真实点估计恢复和 Exact 小预算试跑记录，保存于
`analysis_reproducibility/rare_numerics_20260923`。
