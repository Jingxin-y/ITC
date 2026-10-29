# Indirect treatment comparison: RR / RD / STC

当前版本：2026-10-29 更新，来自 Nibi 的 2026-10-25/26 导出及 10 月 29 日 RD 单条区间修复。

| 目录 | 内容 |
|---|---|
| [code](code) | 与当前结果运行指纹一致的 RR/RD 源码、STC 后处理、非 Wald CI 审阅文件和绘图代码 |
| [results](results) | 最新图表、指标表及核对记录 |
| [original result](original%20result) | 原样逐次 CSV、STC 报告、情景配置和校验信息 |

RR、RD 各有 8 个数据情景 × 3 种分析：`none_MAIC`、`score_unit`、`score_MAIC`；每个分析情景 seeds 1–1000。
8 个数据情景为 common/rare × n=100/500 × 名义 Full ESS 比例 0.15/0.50，固定 diffuse / high / all。n 为每项试验总人数，两个观察臂各 n/2。
RR 的逐次表共 432000 行，RD 共 216000 行。STC RR 单独保存 192000 行后处理结果及配对报告。

- [结果与图表预览](results/README.md)
- [逐次 CSV 检索表](original%20result/config/SCENARIOS.csv)
- [代码入口与对应关系](code/README.md)
- [非 Wald 区间说明](code/CI_methods/README_zh.md)
- [完整性与版本核对](VERIFICATION.json)

`score_unit` 为 score 结局模型使用单位拟合权重；`score_MAIC` 为 score 模型使用 MAIC 拟合权重；`none_MAIC` 不加入 score 并使用 MAIC 拟合权重。三者的目标标准化仍使用 MAIC 权重。
点估计、区间、失败状态及真值均保留。Coverage 以有效返回且有对应真值的区间为分母，另报告区间返回率和占计划重复次数的覆盖比例。
基础逐 seed RDS 仍在集群；本次单条修复的 RDS 保留在 Nibi 原压缩包中，本次操作未同步回集群。此次 GitHub 更新只整理已保存结果和重绘图表，未运行模拟、拟合或区间反演。

2026-10-29 更新：RD `score_unit / S06 / seed 9 / brm-FC-BC` 已采用补算成功的区间；[修复说明](code/repairs/RD_stationary_restart_S06_seed9_20261029/README_zh.md)。RR、RD 仍各为 **3 种设置，每种 8 个参数组合**。旧 none/单位拟合权重补充汇总已移出当前目录。
