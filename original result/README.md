# 最新 Nibi 逐次结果

RR 和 RD 各为 8 个数据情景 × 3 种分析（none_MAIC、score_unit、score_MAIC），每情景1000个seed。
48份 `combined_results.csv` 按导出文件原字节保存，RR每seed18行，RD每seed9行；失败及不适用记录均保留。
RD 的源文件名为 `anchored_rd_estimates.csv`，Nibi 导出时原样命名为 `combined_results.csv`。

- [情景与文件检索](config/SCENARIOS.csv)
- [逐文件字节数和校验值](checksums/FILES.csv)
- [原始导出包SHA256](checksums/SOURCE_ARCHIVES.csv)
- [原始导出配置与来源](provenance)
- [对应图表和指标](../results/README.md)
- [对应代码](../code/README.md)

RR 的 estimate/low/up 是 log(RR)，RD 的 estimate/lower/upper 是自然 RD。
STC RR 放在 [STC_RR](STC_RR)：192000行后处理结果及配对报告，使用已保存的score_unit拟合，未重新拟合。
`combined_results.csv.gz` 与 `paired_results.csv.gz` 是原始大CSV的无损gzip；GitHub普通Git不接受大于100MiB的单文件。校验表同时列出压缩文件和解压后原CSV的SHA256。
STC的两个小RDS仅保存配置/运行元信息，不包含逐seed模型对象或个体数据；其余原始RDS仍在集群。

provenance 中的 exported_file 字段沿用原压缩包内路径。当前仓库导航使用 config/SCENARIOS.csv 的 csv_path。
provenance 内原 MANIFEST.sha256 用于核对原导出包；当前目录请使用本目录 SHA256SUMS 和 checksums/FILES.csv。
完整性核对见根目录 [VERIFICATION.json](../VERIFICATION.json)。
