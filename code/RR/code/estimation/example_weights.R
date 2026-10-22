## 在 anchored_v1 工作目录运行：先生成数据，再选择协变量求 MAIC 权重。
source("data/example_anchored.R", encoding = "UTF-8")
source("estimation/helpers/common.R", encoding = "UTF-8")
source("estimation/anchored_weights.R", encoding = "UTF-8")

## 在调用 get.weight() 前选择协变量，组成矩阵和对应的目标均值向量。
weight_covariates <- c("z1", "z2", "z3")  # 示例：只平衡真正的 effect modifiers
# weight_covariates <- paste0("z", 1:6)   # 平衡全部六个变量
# weight_covariates <- c("z1", "z4")      # 任意组合
# weight_covariates <- character(0)      # 不平衡协变量，全部权重为 1
Z.weight <- as.matrix(data.IPD$data[, weight_covariates, drop = FALSE])
target.mean <- data.AD$anchored.mean.x[weight_covariates]
weight.fit <- get.weight(Z.weight, target.mean)
weight.fit$covariates
weight.fit$ga
## 仍展示六个均值；只要求所选协变量的加权均值匹配目标。
rbind(weighted_AC = colSums(covariate_matrix(data.IPD$data) * weight.fit$weight) /
                      sum(weight.fit$weight),
      sample_BC = data.AD$anchored.mean.x,
      population_BC = design$e)
sum(weight.fit$weight)^2 / sum(weight.fit$weight^2) / nrow(data.IPD$data)

## 若选择已经生成的分数列，先取得对应的 AD 汇总，再按这些列求权重：
# score.AD <- get.aggre(data.AD.or, covariates = "h_E")
# Z.score <- as.matrix(data.IPD$data[, "h_E", drop = FALSE])
# score.weight <- get.weight(Z.score, score.AD$anchored.mean.x)

## 后续联合 sandwich 请直接传 weight.fit$ga（保留变量名称），
## 或显式指定 weight_covariates=weight.fit$covariates；结局模型另行设置。
