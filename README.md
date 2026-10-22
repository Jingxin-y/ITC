# Indirect treatment comparison: RR / RD

RR 与 RD 的源码、当前 result_v2 绘图结果和逐次估计 CSV。

| 目录 | 内容 |
|---|---|
| [code](code) | 完整 RR/RD 源码、非 Wald CI 审阅文件、绘图程序和配置 |
| [results](results) | RR、RD 各 288 个情景的 PNG、PDF、指标表和说明 |
| [original result](original%20result) | RR、RD 各 288 份逐次结果 CSV、情景配置和校验信息 |

- [逐次原始结果与校验](original%20result/README.md)
- [中文代码使用说明](code/README_zh.txt)
- [非 Wald 置信区间](code/CI_methods/README_zh.md)
- [绘图代码](code/plotting/README_zh.txt)
- [RR 结果](results/RR) · [RD 结果](results/RD) · [情景检索表](results/SCENARIOS.csv)

代码对应 `anchored_v1/result_cluster/result_v2` 的模型和算法；结果文件复制自当前 `plots_RR`、`plots_RD`。
每尺度有 288 个情景，每情景 MC1000。逐 seed 的完整 RDS 未包含在此仓库。

## Quick code checks

```sh
cd code
python verify_checksums.py
Rscript --vanilla verify_runtime.R RR
Rscript --vanilla verify_runtime.R RD
```

One-seed interface examples and plotting commands are documented inside `code/`.
Existing upstream license notices are preserved; no new project-wide license is assigned here.
