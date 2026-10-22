## Calculate MAIC weights first; fit models with unit weights by default.
## Marginal standardization and CMH stratum counts use the MAIC weights.
get.compare.anchored <- function(data.ipd, data.ad,
                                 param = "RR",
                                 run_exact = TRUE,
                                 seed = NULL, exact_seed_offset = 1000000L,
                                 bc_control = list(), weight_covariates = NULL,
                                 model_covariates = NULL, effect_modifiers = NULL,
                                 return_full = length(effect_modifiers) > 0L,
                                 non_effect_modifiers = NULL,
                                 op_model = c("marginal", "linear"), score_spec = NULL,
                                 score_standardization = c("bc_estimated", "fixed"),
                                 fit_control = list(), bc_target = c("marginal", "coefficient"),
                                 model_weighting = c("unit", "maic")) {
  model_weighting <- match.arg(model_weighting)
  use_maic_in_model <- model_weighting == "maic"
  bc_target <- match.arg(bc_target)
  if (bc_target == "marginal") validate_marginal_bc_control(bc_control)
  fit_control <- anchored_fit_control(fit_control)
  if (!identical(param, "RR")) stop("The comparison supports RR only; use get.estimate() for an individual RD fit.")
  score_standardization <- match.arg(score_standardization)
  if (is.null(score_spec)) score_spec <- data.ipd$score_spec
  if (is.null(score_spec)) score_spec <- attr(data.ipd$data,"brm_score_spec")
  uses_score <- any(c(model_covariates,effect_modifiers,non_effect_modifiers,weight_covariates) %in%
                    c("S_E","S_P","h_E","h_P"))
  if (is.null(score_spec) && uses_score && score_standardization == "bc_estimated") {
    score_spec <- make_brm_score_spec()
    if (is.null(effect_modifiers)) effect_modifiers <- character(0)
    if (is.null(model_covariates)) model_covariates <- covariate_names(data.ipd$data)
    if (is.null(non_effect_modifiers)) non_effect_modifiers <- setdiff(model_covariates,effect_modifiers)
  }
  if (!is.null(score_spec)) {
    if (score_standardization == "bc_estimated") score_spec <- brm_scores_from_bc(score_spec,data.ad)
  }
  prepared <- prepare_brm_score_model(data.ipd, score_spec, effect_modifiers,
                                      non_effect_modifiers, model_covariates)
  data.ipd <- prepared$data
  score_spec <- prepared$score_spec
  effect_modifiers <- prepared$effect_modifiers
  non_effect_modifiers <- prepared$non_effect_modifiers
  model_covariates <- prepared$model_covariates
  if (!is.null(score_spec)) {
    score_weights <- intersect(weight_covariates,c("S_E","S_P","h_E","h_P"))
    data.ad <- brm_score_summary(data.ad,score_spec,score_weights)
  }
  y <- data.ipd$data$y
  if (is.null(model_covariates)) model_covariates <- covariate_names(data.ipd$data)
  op_model <- match.arg(op_model)
  ## 输入错误在 safe_try() 外报告，不能被转换成整行 NA。
  op.design <- brm_op_design(data.ipd$data, effect_modifiers, non_effect_modifiers,
                            model_covariates, op_model)
  op.spec <- op.design$spec
  non_effect_modifiers <- op.spec$non_effect_modifiers
  if (length(effect_modifiers)) return_full <- TRUE # 多参数模型始终保留完整结果。

  Z.weight <- covariate_matrix(data.ipd$data, weight_covariates)
  target.mean.weight <- select_weight_summary(data.ad$anchored.mean.x,
    colnames(Z.weight), covariate_names(data.ipd$data))
  model.anweight <- tryCatch(get.weight(Z.weight, target.mean.weight),
    maic_calibration_error = function(e) {
      warning(conditionMessage(e), call. = FALSE)
      e$result
    })
  model.weights <- if (use_maic_in_model) model.anweight else list(
    weight = rep(1, length(y)), ga = setNames(numeric(0), character(0)),
    covariates = character(0), target_mean = numeric(0), status = "available",
    diagnostics = list(mode = "fixed_unit", estimated_weight_uncertainty = FALSE))
  anchored.weight <- model.weights$weight
  anchored.gamma <- model.weights$ga
  model.weight.covariates <- model.weights$covariates
  ## "AD" selects the existing fixed-weight covariance branch; the observations
  ## passed to get.estimate() remain the same AC individual data.
  model.group <- if (use_maic_in_model) "IPD" else "AD"

  variables <- list(weight_covariates = model.anweight$covariates,
    model_covariates = model_covariates, effect_modifiers = effect_modifiers,
    non_effect_modifiers = non_effect_modifiers, op_model = op_model, score_spec = score_spec,
    score_standardization = score_standardization, fit_control = fit_control, bc_target = bc_target,
    model_weighting = model_weighting, standardization_weighting = "maic",
    cmh_weighting = "maic",
    cmh_sandwich_uncertainty = if (use_maic_in_model) "joint" else "fixed",
    comparison_type = "anchored", method_naming_version = 2L)
  schema <- comparison_schema(data.ipd$data, model_covariates, effect_modifiers, op.design)
  if (model.weights$status != "available") {
    fits <- setNames(vector("list", length(anchored_method_names())), anchored_method_names())
    intervals <- setNames(rep(list(make_na_exact(ncol(op.design$va))), 3),
                          c("brm.bc", "firth.bc", "bayes.bc"))
    if (bc_target == "marginal") intervals <- lapply(intervals,function(x)
      make_na_marginal_exact("weighting_failed"))
    return(finish_anchored_comparison(fits, intervals, model.anweight, variables, schema,
                                     run_exact, return_full, model_weights = model.weights))
  }
  anchored.est <- safe_try(
    get.estimate(model.group, data.ipd, anchored.weight, anchored.gamma,
                 target_mean = data.ad$anchored.mean.x,
                 target_var = data.ad$anchored.var.x,
                 m_AD = data.ad$anchored.n,
                 target_cov = data.ad$anchored.cov.x,
                 weight_covariates = model.weight.covariates,
                 model_covariates = model_covariates, effect_modifiers = effect_modifiers,
                 non_effect_modifiers = non_effect_modifiers, op_model = op_model,
                 fit_control = fit_control),
    "brm"
  )
  if (is.null(anchored.est)) anchored.est <- make_na_est()

  LB.anchored <- safe_try(fit_logbinomial(data.ipd$data, anchored.weight,
                                          "LB.anchored", model_covariates,
                                          effect_modifiers, maxit = fit_control$glm_maxit), "LB")
  RLB.sandwich.anchored <- safe_try(
    fit_robust_lb_joint(data.ipd$data, anchored.weight,
      gamma = anchored.gamma,
      target_mean = data.ad$anchored.mean.x,
      target_var = data.ad$anchored.var.x,
      m_AD = data.ad$anchored.n,
      target_cov = data.ad$anchored.cov.x,
      account_for_gamma = use_maic_in_model,
      account_for_target_mean = use_maic_in_model,
      weight_covariates = model.weight.covariates,
      model_covariates = model_covariates, effect_modifiers = effect_modifiers,
      return_full = TRUE, maxit = fit_control$glm_maxit), "RLB")
  LP.anchored <- safe_try(fit_lp(data.ipd$data, anchored.weight,
                                 model_covariates, effect_modifiers, maxit = fit_control$glm_maxit), "LP")
  RLP.anchored <- safe_try(
    fit_robust_lp_joint(data.ipd$data, anchored.weight,
                        gamma = anchored.gamma,
                        target_mean = data.ad$anchored.mean.x,
                        target_var = data.ad$anchored.var.x,
                        m_AD = data.ad$anchored.n,
                        target_cov = data.ad$anchored.cov.x,
                        account_for_gamma = use_maic_in_model,
                        account_for_target_mean = use_maic_in_model,
                        weight_covariates = model.weight.covariates,
                        model_covariates = model_covariates, effect_modifiers = effect_modifiers,
                        return_full = TRUE, maxit = fit_control$glm_maxit),
    "RLP"
  )
  ## CMH has no separate regression fit/standardization step: transport its
  ## stratum counts even when regression outcome models use unit weights.
  ## Keep the selected fixed/joint uncertainty scope of the comparison.
  CMH.anchored <- CMH.sandwich.anchored <- NULL
  if (model.anweight$status == "available") {
    CMH.anchored <- safe_try(get.weighted.RR(data.ipd$data, model.anweight$weight,
                                             model_covariates = model_covariates),
                             "CMH")
    CMH.sandwich.anchored <- safe_try(
      get.weighted.RR.sandwich(data.ipd$data, model.anweight$weight,
                             gamma = model.anweight$ga,
                             target_mean = data.ad$anchored.mean.x,
                             target_var = data.ad$anchored.var.x,
                             m_AD = data.ad$anchored.n,
                             target_cov = data.ad$anchored.cov.x,
                             account_for_gamma = use_maic_in_model,
                             account_for_target_mean = use_maic_in_model,
                             weight_covariates = model.anweight$covariates,
                             model_covariates = model_covariates),
      "CMH_sandwich"
    )
  }

  bayes.anchored <- if (length(effect_modifiers)) NULL else
    safe_try(get.brm.ad(anchored.est, data.ipd$data, anchored.weight), "brmad")
  if (is.null(bayes.anchored)) bayes.anchored <- make_na_est()

  glm.fc <- fit_anchored_glm_fc(data.ipd$data, anchored.weight,
    model_covariates, effect_modifiers, fit_control, use_maic_in_model,
    joint_args = list(gamma = anchored.gamma, target_mean = data.ad$anchored.mean.x,
      target_var = data.ad$anchored.var.x, m_AD = data.ad$anchored.n,
      target_cov = data.ad$anchored.cov.x, weight_covariates = model.weight.covariates))

  brm.bc <- get.brm.bc(
    data.ipd, data.ad, anchored.weight, anchored.gamma,
    anchored.est, bayes.anchored, run_exact = run_exact,
    seed = seed, exact_seed_offset = exact_seed_offset, bc_control = bc_control,
    weight_covariates = model.weight.covariates,
    model_covariates = model_covariates, effect_modifiers = effect_modifiers,
    non_effect_modifiers = non_effect_modifiers, op_model = op_model, fit_control = fit_control,
    bc_target = bc_target, model_weighting = model_weighting,
    standardization.weight = model.anweight$weight
  )

  fits <- list(brm = anchored.est, LB = LB.anchored,
    LP = LP.anchored, RLP = RLP.anchored, CMH = CMH.anchored,
    CMH_sandwich = CMH.sandwich.anchored, brmad = bayes.anchored,
    "brm-BC" = anchored.est, "brm-FC-BC" = brm.bc$firth, "brm_b-BC" = bayes.anchored,
    "brm-FC" = brm.bc$firth, RLB = RLB.sandwich.anchored)
  fits <- c(fits, glm.fc)
  finish_anchored_comparison(fits, brm.bc[c("brm.bc", "firth.bc", "bayes.bc")],
                             model.anweight, variables, schema, run_exact, return_full,
                             model_weights = model.weights)
}
run.anchored <- function(design, n_AC, n_BC = n_AC,
                         run_exact = FALSE, seed = NULL,
                         exact_seed_offset = 1000000L, bc_control = list(),
                         weight_covariates = NULL, model_covariates = NULL,
                         effect_modifiers = NULL,
                         return_full = TRUE,
                         non_effect_modifiers = NULL,
                         op_model = c("marginal", "linear"), score_spec = NULL,
                         return_data = FALSE,
                         score_standardization = c("bc_estimated", "fixed"),
                         fit_control = list(), bc_target = c("marginal", "coefficient"),
                         model_weighting = c("unit", "maic")) {
  model_weighting <- match.arg(model_weighting)
  bc_target <- match.arg(bc_target)
  if (!is.list(design) || !all(c("parameters", "e", "q", "risk_table", "marginal_truth") %in% names(design)))
    stop("design must be returned by webster_rr_design() or joint_rr_gop_design(); the legacy event/hypothesis API is retired.")
  for (size in list(n_AC, n_BC)) {
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size < 2 || size %% 2 != 0)
      stop("n_AC and n_BC must be positive even integers.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max - exact_seed_offset)
      stop("seed must be a nonnegative integer compatible with exact_seed_offset.")
    set.seed(seed)
  }
  sim <- generate_anchored_data(n_AC, n_BC, design)
  ## 无自定义配置时也允许按现成的分数列平衡；自定义分数由比较入口转换 AD 原始矩。
  if (is.null(score_spec) && length(weight_covariates)) {
    summary_cols <- union(covariate_names(sim$data.AD.or$data), weight_covariates)
    sim$data.AD <- get.aggre(sim$data.AD.or, covariates = summary_cols)
  }
  result <- get.compare.anchored(
    sim$data.IPD, sim$data.AD,
    param = "RR",
    run_exact = run_exact,
    seed = seed, exact_seed_offset = exact_seed_offset, bc_control = bc_control,
    weight_covariates = weight_covariates, model_covariates = model_covariates,
    effect_modifiers = effect_modifiers, return_full = return_full,
    non_effect_modifiers = non_effect_modifiers, op_model = match.arg(op_model),
    score_spec = score_spec, score_standardization = match.arg(score_standardization),
    fit_control = fit_control, bc_target = bc_target, model_weighting = model_weighting
  )
  metadata <- list(design = design, n_AC = n_AC, n_BC = n_BC, seed = seed)
  if (is.list(result)) {
    result$simulation <- metadata
    result$reference_truth <- brm_score_coefficient_truth(design,result$variables)
    if (return_data) result$data <- sim[c("data.IPD", "data.AD.or", "data.AD")]
  } else {
    attr(result, "simulation") <- metadata
    if (return_data) attr(result, "data") <- sim[c("data.IPD", "data.AD.or", "data.AD")]
  }
  result
}
