## 演示配置；在 anchored_v1 运行：
## Rscript estimation/run_simulation_anchored.R estimation/simulation_config.R results_unit_fit_maic_target
## 数据真值与分析假设分开。修改 scores 不会改变生成数据。
## Webster-Clark RR log-risk adaptation; still z1:z3 EM and z4:z6 non-EM.
## Defaults recover original AC and BC event rates at ESS=.5; A/C differs from B/C.
## RD generation is available separately; this configuration remains an RR analysis.
design <- webster_rr_design(event = "common", ess_ratio = .5)
scores <- make_brm_score_spec(
  em_covariates = c("z1", "z2", "z3"),
  non_em_covariates = c("z4", "z5", "z6"))
simulation_args <- list(
  design = design, n_AC = 1000, n_BC = 1000, seeds = 1:10,
  score_spec = scores, weight_covariates = paste0("z", 1:6),
  model_weighting = "unit", # Fit all models with 1; MAIC weights still standardize marginal effects.
  op_model = "marginal", run_exact = FALSE, bc_target = "marginal", score_standardization = "bc_estimated",
  fit_control = list(brm_max_step=1000L,brm_optim_maxit=200L,fc_max_step=300L,glm_maxit=200L),
  bc_control = list(optim.maxit=200L))

## run_exact=TRUE 显式调用 exact_marginal_target()；输出另存 anchored_joint_marginal_exact.csv。
## marginal 行的覆盖率需要显式提供 marginal logRR 真值，不能使用条件 alpha 真值。

## NULL: 对能准确表示真实 RR 的 BC 标准化分数，自动使用每次重复对应的系数真值。
## 每次 BC 的样本中心和尺度不同，不能机械沿用固定 alpha/lambda。
## 自定义真值可用 r/method/component/term/truth 五列；r 与 seeds 一致。
## 错误 EM 分数若不能准确表示真实 RR，不自动给条件系数填写真值。
truth <- NULL
