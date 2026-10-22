# Original per-replicate results: RR / RD

与当前 `result_v2` 绘图对应的逐次估计 CSV。RR、RD 各 288 个情景，每情景保留 seeds 1–1000；包含失败记录及其状态。

| 目录 | 文件与用途 |
|---|---|
| [RR](RR) | 每情景 `combined_results.csv`：每 seed 18 条方法/区间记录 |
| [RD](RD) | 每情景 `anchored_rd_estimates.csv`：每 seed 9 条方法/区间记录 |
| [config/SCENARIOS.csv](config/SCENARIOS.csv) | 情景参数、seed 数量、逐次 CSV 相对路径和集群来源 |
| [checksums/FILES.csv](checksums/FILES.csv) | 全部 576 份 CSV 的字节数、MD5、SHA256 与来源 |
| [SHA256SUMS](SHA256SUMS) | 原精简包内文件的 SHA256 |
| [VERIFICATION.json](VERIFICATION.json) | 上传前的独立完整性核对记录 |

RR 的 `estimate/low/up` 使用 log(RR) 尺度；`Exact` 是 BC 区间的原导出标签。RD 的 `estimate/lower/upper` 使用自然 RD 尺度。
CSV 保留真值、点估计、区间和有效状态，供核对 coverage、Bias、SE、Accuracy 及重新绘图。
Coverage 的分母为有效返回且有对应真值的区间数，需同时查看区间返回比例。

来源为集群精简打包任务 `24971435`，沿用原始任务 `24947571` 的冻结清单。全部 CSV 按原字节保存；校验了各情景全部 1000 个 seed、方法/区间记录、原清单 MD5 和文件 SHA256。
文件中的 scratch 路径用于记录来源；本目录的导航使用相对路径。完整 RDS、模型对象和模拟个体数据仍保存在集群。
情景 ESS 比较允许 `1e-9` 的数值尾差；原 CSV 和配置文件均未改写。

- [对应绘图结果](../results)
- [计算及区间代码](../code)
- [绘图运行说明](../code/plotting/README_zh.txt)

`SHA256SUMS` 是原精简包的清单，不覆盖后加的本说明、`.gitattributes` 和 `VERIFICATION.json`。
