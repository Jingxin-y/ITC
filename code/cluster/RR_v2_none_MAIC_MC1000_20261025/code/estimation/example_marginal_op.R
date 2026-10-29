## 在 anchored_v1 根目录逐段运行：用户指定 EM 和非 EM。
source("data/example_anchored.R", encoding="UTF-8")
source_anchored()

## 在估计阶段自选原始变量；真实生成机制不变。
scores <- make_brm_score_spec(
  em_covariates=c("z1","z2","z3"),
  non_em_covariates=c("z4","z5","z6")
)
em <- "h_E"
non_em <- "h_P"
weight_cols <- paste0("z", 1:6)
wf <- get.weight(covariate_matrix(data.IPD$data, weight_cols),
                 data.AD$anchored.mean.x[weight_cols])

fit <- get.estimate(
  "IPD", data.IPD, wf$weight, wf$ga,
  target_mean=data.AD$anchored.mean.x,
  target_cov=data.AD$anchored.cov.x, m_AD=data.AD$anchored.n,
  effect_modifiers=em, non_effect_modifiers=non_em,
  op_model="marginal", score_spec=scores
)
fit$coefficients
fit$op_spec$cells
## log RR(h_E)=alpha[1]+alpha[2]*h_E；OP 对 h_E 的观察值饱和。
## 不使用 h_P 的回归项，也不使用生成器的 risk_table 进行估计。

## 完全 marginal：两部分均只有截距，目标为当前加权人群的 marginal RR。
fully_marginal <- get.estimate(
  "IPD", data.IPD, wf$weight, wf$ga,
  target_mean=data.AD$anchored.mean.x,
  target_cov=data.AD$anchored.cov.x, m_AD=data.AD$anchored.n,
  effect_modifiers=character(0), non_effect_modifiers=c(em, non_em), score_spec=scores
)
fully_marginal$coefficients

## 保留原来的线性 OP 作为错设/敏感性对照。
linear_op <- get.estimate(
  "IPD", data.IPD, wf$weight, wf$ga,
  target_mean=data.AD$anchored.mean.x,
  target_cov=data.AD$anchored.cov.x, m_AD=data.AD$anchored.n,
  effect_modifiers=em, non_effect_modifiers=non_em,
  model_covariates=c(em, non_em), op_model="linear", score_spec=scores
)
