## 将模拟的 AD 个体数据转成协变量汇总信息。
## 原代码：compare_RR_anchored.R

get.aggre <- function(data, covariates = NULL) {
  Z <- covariate_matrix(data$data, covariates)
  m.AD <- nrow(Z)
  q <- ncol(Z)

  ## 与按两臂样本量加权的均值相同，也适用于某臂没有受试者时。
  anchored.mean <- colMeans(Z)
  anchored.cov <- if (m.AD > 1) {
    as.matrix(stats::cov(Z))
  } else {
    matrix(NA_real_, nrow = q, ncol = q,
           dimnames = list(colnames(Z), colnames(Z)))
  }
  anchored.var <- diag(anchored.cov)

  list(
    anchored.mean.x = anchored.mean,
    anchored.var.x = anchored.var,
    anchored.cov.x = anchored.cov,
    anchored.n = m.AD
  )
}

