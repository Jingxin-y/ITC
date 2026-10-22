## 加载三 EM / 三非 EM 的直接 RR、独立 RD 及联合 RR-GOP 生成机制。
source_anchored_data <- function(code_dir = ".") {
  files <- c(
    "data/helpers/covariates.R",
    "data/helpers/getProbRR_GOP.R",
    "data/helpers/joint_rr_gop.R",
    "data/helpers/webster_design.R",
    "data/ad_summary.R",
    "data/generate_RR.R",
    "data/generate_webster_RR.R",
    "data/generate_webster_RR_scenarios.R",
    "data/generate_webster_RD.R"
  )
  for (file in files) source(file.path(code_dir, file), encoding = "UTF-8")
  invisible(NULL)
}

## 加载当前联合 RR-GOP 生成、估计与批量接口。
source_anchored <- function(code_dir = ".",
                            exact_file = "CI_exact_adapt_para.R", legacy = FALSE) {
  if (isTRUE(legacy)) stop("The legacy generator is no longer included. Use joint_rr_gop_design() and run.anchored(design, n_AC, n_BC).")
  source_anchored_data(code_dir)
  library(MASS)

  files <- c(
    "estimation/helpers/common.R",
    "estimation/helpers/brm_scores.R",
    "estimation/helpers/brm_op_design.R",
    "estimation/helpers/comparison_status.R",
    "estimation/helpers/getProbScalarRR.R",
    "estimation/helpers/getProbScalarRD.R",
    "estimation/helpers/1.1_MLE_Point.R",
    "estimation/helpers/1.2_MLE_Var.R",
    "estimation/helpers/1.2_MLE_Var_joint.R",
    "estimation/helpers/firth_augmentation.R",
    "estimation/helpers/MLE_Point_Firth.R",
    "estimation/helpers/MLE_Point_Firth_for_RR.R",
    "estimation/helpers/MLE_Point_Firth_for_RD.R",
    "estimation/helpers/1_CallMLE.R",
    "estimation/anchored_weights.R",
    "estimation/anchored_regression.R",
    "estimation/anchored_glm_fc.R",
    "estimation/anchored_cmh.R",
    "estimation/anchored_bayes_gcomp.R",
    "estimation/anchored_comparison.R",
    "estimation/anchored_simulation.R",
    "summary/summarize_joint_simulation.R"
  )
  for (file in files) source(file.path(code_dir, file), encoding = "UTF-8")

  source(file.path(code_dir, "estimation", "helpers", exact_file), encoding = "UTF-8")
  source(file.path(code_dir, "estimation/helpers/CI_linear_target.R"), encoding = "UTF-8")
  source(file.path(code_dir, "estimation/helpers/CI_marginal_target.R"), encoding = "UTF-8")
  source(file.path(code_dir, "estimation/helpers/marginal_comparison.R"), encoding = "UTF-8")

  invisible(NULL)
}
