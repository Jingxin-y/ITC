## 当前 DGM 的重复模拟；共享 design，保留完整参数及每次拟合状态。
anchored_method_names <- function() c(
  "brm", "LB", "LP", "RLP", "CMH", "CMH_sandwich",
  "brmad", "brm-BC", "brm-FC-BC", "brm_b-BC", "brm-FC", "RLB",
  anchored_glm_fc_method_names())

simulate_anchored <- function(design, n_AC, n_BC = n_AC, seeds = seq_len(100),
                              run_exact = FALSE, weight_covariates = NULL,
                              model_covariates = NULL, effect_modifiers = NULL,
                              non_effect_modifiers = NULL, score_spec = NULL,
                              op_model = c("marginal", "linear"),
                              exact_seed_offset = 1000000L, bc_control = list(),
                              checkpoint_file = NULL,
                              score_standardization = c("bc_estimated", "fixed"),
                              fit_control = list(), bc_target = c("marginal", "coefficient"),
                              model_weighting = c("unit", "maic")) {
  model_weighting <- match.arg(model_weighting)
  bc_target <- match.arg(bc_target)
  if (!is.numeric(seeds) || !length(seeds) || any(!is.finite(seeds)) ||
      any(seeds < 0 | seeds != floor(seeds)) || anyDuplicated(seeds))
    stop("seeds must be distinct nonnegative integers.")
  op_model <- match.arg(op_model)
  score_standardization <- match.arg(score_standardization)
  fit_control <- anchored_fit_control(fit_control)
  rows <- diagnostics <- vector("list", length(seeds))
  reference_truth <- score_references <- vector("list",length(seeds))
  nonconverged_rows <- nonconverged_fits <- vector("list",length(seeds))
  variables <- NULL
  snapshot <- function(completed) structure(list(
    estimates = do.call(rbind, rows[seq_len(completed)]), diagnostics = diagnostics[seq_len(completed)],
    nonconverged_estimates = do.call(rbind,nonconverged_rows[seq_len(completed)]),
    nonconverged_fits = Filter(Negate(is.null),nonconverged_fits[seq_len(completed)]),
    parameter_schema = if (completed) unique(do.call(rbind, rows[seq_len(completed)])[
      c("method", "component", "parameter", "term")]) else NULL,
    design = design, n_AC = n_AC, n_BC = n_BC, seeds = seeds[seq_len(completed)],
    requested_seeds = seeds, complete = completed == length(seeds),
    reference_truth = do.call(rbind,reference_truth[seq_len(completed)]),
    score_references = score_references[seq_len(completed)],
    score_standardization = score_standardization, model_weighting = model_weighting,
    run_exact = run_exact, variables = variables, fit_control = fit_control, bc_target = bc_target,
    exact_seed_offset = exact_seed_offset, bc_control = bc_control),
    class = c("anchored_simulation", "list"))
  if (!is.null(checkpoint_file)) {
    dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(snapshot(0L), checkpoint_file)
  }
  for (i in seq_along(seeds)) {
    messages <- character(0)
    fit <- tryCatch(withCallingHandlers(run.anchored(
      design = design, n_AC = n_AC, n_BC = n_BC, seed = seeds[i], run_exact = run_exact,
      weight_covariates = weight_covariates, model_covariates = model_covariates,
      effect_modifiers = effect_modifiers, non_effect_modifiers = non_effect_modifiers,
      score_spec = score_spec, op_model = op_model, return_full = TRUE,
      score_standardization = score_standardization, fit_control = fit_control,
      exact_seed_offset = exact_seed_offset, bc_control = bc_control, bc_target = bc_target,
      model_weighting = model_weighting),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        if (!is.null(checkpoint_file)) {
          partial <- snapshot(i - 1L)
          partial$interrupted_at_seed <- seeds[i]
          partial$error <- conditionMessage(e)
          saveRDS(partial, checkpoint_file)
        }
        stop(e)
      })
    rows[[i]] <- cbind(r = seeds[i], fit$estimates)
    raw <- fit$nonconverged_estimates
    nonconverged_rows[[i]] <- cbind(r=rep(seeds[i],nrow(raw)),raw)
    if (length(fit$nonconverged_fits)) nonconverged_fits[[i]] <- list(
      r=seeds[i],fits=fit$nonconverged_fits,weights=fit$weights,
      model_weights=fit$model_weights,variables=fit$variables)
    if (!is.null(fit$reference_truth)) reference_truth[[i]] <- cbind(r=seeds[i],fit$reference_truth)
    score_references[[i]] <- list(r=seeds[i],score_spec=fit$variables$score_spec)
    w <- fit$weights$weight
    mw <- fit$model_weights$weight
    diagnostics[[i]] <- list(r = seeds[i], warnings = unique(messages),
      ess = if (all(is.finite(w)) && sum(w^2) > 0) sum(w)^2 / sum(w^2) else NA_real_,
      model_ess = if (all(is.finite(mw)) && sum(mw^2) > 0) sum(mw)^2 / sum(mw^2) else NA_real_,
      model_weighting = model_weighting, model_weight_status = fit$model_weights$status,
      weights = fit$weights$diagnostics, weight_status = fit$weights$status,
      fit_status = unique(fit$estimates[c("method", "status")]),
      brm = lapply(fit$fits[c("brm", "brm-FC")], function(f) {
        list(converged = f$converged, op_spec = f$op_spec, fc = f$fc.diagnostics,
             optimization = f$optimization)
      }), glm_fc = lapply(fit$fits[anchored_glm_fc_method_names()], function(f)
        list(converged = f$converged, fc = f$fc.diagnostics, variance = f$variance.method)),
      bc = lapply(fit$bc, function(x) x$diagnostics))
    variables <- fit$variables
    variables$score_spec <- score_spec # selection template; realized units are in score_references
    if (!is.null(checkpoint_file)) saveRDS(snapshot(i), checkpoint_file)
  }
  snapshot(length(seeds))
}
