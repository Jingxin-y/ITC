# 当前运行包快照

四个包均按原 `SHA256SUMS` 所列文件原样保存，包含自身 `code/`、配置和运行入口。

- `RR_v2_full24_MC1000_20261024`：RR 的 score_unit / score_MAIC 来源；原 full24 还包含 none（单位拟合权重）。
- `RR_v2_none_MAIC_MC1000_20261025`：当前 RR none_MAIC 来源。
- `RD_v2_full24_MC1000_20261024`：RD 的 score_unit / score_MAIC 来源；原 full24 还包含 none（单位拟合权重）。
- `RD_v2_noneMAIC_MC1000_20261025`：当前 RD none_MAIC 来源。

包内 README 和验证说明保留准备时的原文与状态记录。当前已导出结果的完成情况和代码指纹核对，以仓库根目录 `VERIFICATION.json` 为准。
当前主结果为每尺度 8 个数据情景 × none_MAIC / score_unit / score_MAIC 三种分析。
提交会重新运行模拟，请仅在确实需要重跑时，按包内说明上传到集群并提交。
