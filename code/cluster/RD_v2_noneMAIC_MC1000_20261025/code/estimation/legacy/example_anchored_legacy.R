## 旧数据生成机制下的 11 方法比较示例。
## 先将 R 的工作目录设为 anchored_v1，然后逐段运行。
source("load_anchored.R", encoding = "UTF-8")
source_anchored(legacy = TRUE)

## 1. 生成 A-C IPD 和 B-C AD。
set.seed(1)
data.IPD <- data.generation_legacy(
  "RR", "IPD", n = 100, event = "common", hypothesis = "alternative",
  n_cov = 3, ipd_prob = c(0.7, 0.3, 0.7)
)
data.AD.or <- data.generation_legacy(
  "RR", "AD", n = 100, event = "common", hypothesis = "alternative",
  n_cov = 3, ipd_prob = c(0.7, 0.3, 0.7),
  ess_ratio = 0.5, direction = c(-1, 1, -1)
)

## 2. 提取 AD 的均值、方差、协方差和样本量。
data.AD <- get.aggre(data.AD.or)

## 3. 计算 MAIC 权重，可直接查看加权后的协变量均值。
weight.fit <- get.weight(covariate_matrix(data.IPD$data), data.AD$anchored.mean.x)
colSums(covariate_matrix(data.IPD$data) * weight.fit$weight) / sum(weight.fit$weight)
data.AD$anchored.mean.x

## 4. Anchored 比较（包括三种 BRM-BC）。BC 较慢，此阅读示例先关闭；正式运行可改为 TRUE。
result <- get.compare.anchored(
  data.IPD, data.AD,
  run_exact = FALSE, seed = 1
)
print(result)
