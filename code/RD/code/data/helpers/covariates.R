## 提取原始协变量；数据汇总与估计共用，先由 source_anchored_data() 加载。
covariate_names <- function(data) {
  grep("^z[0-9]+$", names(data), value = TRUE)
}

covariate_matrix <- function(data, covariates = NULL) {
  z.cols <- if (is.null(covariates)) covariate_names(data) else covariates
  as.matrix(data[, z.cols, drop = FALSE])
}
