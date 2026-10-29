rd_empty_rows <- function(model_weighting = "unit") {
  out <- rd_method_plan(); out$component <- "marginal"; out$parameter <- "RD"
  out$term <- "RD_AC_in_BC"; out$effect_scale <- "RD"; out$interval <- out$ci_type
  out$target <- "BC population marginal A-C risk difference"
  for (key in c("estimate", "se", "lower", "upper", "p.value", "risk_C", "risk_A",
                "variance_ipd", "variance_target", "n_eff_C", "n_eff_A", "se_from_ci_width"))
    out[[key]] <- NA_real_
  out$status <- out$se_status <- out$interval_status <- out$p_status <- "pending"
  out$message <- out$approximation <- ""; out$counterfactual_risks_valid <- NA
  out$fitting_weights <- ifelse(out$method %in% c("MN_eff", "RD_MAIC_sandwich") |
    model_weighting == "maic", "MAIC", "ones")
  out$standardization_weights <- "MAIC"
  out$uncertainty <- if (model_weighting == "maic")
    "joint weights, empirical standardization, AD means" else "fixed fitting/standardization weights and covariates"
  out$uncertainty[out$ci_type == "BC"] <- "fixed covariates, treatment and both weight sets; parametric bootstrap"
  out$uncertainty[out$method == "MN_eff"] <- "fixed MAIC weights; Kish effective counts approximation"
  out
}
validate_rd_bc_control <- function(control) {
  reserved <- c("param", "y", "x", "va", "vb", "weight", "pars", "pa", "pb",
    "va.target", "vb.target", "standardization.weight", "vcov.pars", "psi.se", "seed")
  allowed <- setdiff(names(formals(exact_marginal_rd_target)), reserved)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% allowed))))
    stop("Unknown, duplicate or reserved RD bc_control arguments.")
}
rd_brm_basis <- function(ipd, score_spec,
    brm_effect_model = c("linear_atanh", "saturated")) {
  brm_effect_model <- match.arg(brm_effect_model)
  selected <- score_spec$em_covariates
  if (!length(selected)) levels <- 0 else {
    combos <- expand.grid(rep(list(0:1), length(selected)))
    total <- rep(0, nrow(combos))
    for (j in seq_along(selected)) {
      z <- combos[[j]]; if (selected[j] %in% score_spec$reverse_covariates) z <- 1 - z
      total <- total + score_spec$coefficients[selected[j]] * z
    }
    levels <- sort(unique(total))
  }
  if (length(levels) <= 1L) {
    cols <- character(0); model <- "intercept_only"; label <- "intercept only"
  } else if (brm_effect_model == "linear_atanh") {
    if (!"h_E" %in% names(ipd$data) || any(!is.finite(ipd$data$h_E)))
      stop("Linear atanh-RD requires a finite prepared h_E score.")
    cols <- "h_E"; model <- "linear_atanh"; label <- "linear h_E on atanh(RD)"
  } else {
    cols <- paste0("RD_EM_", seq.int(2L, length(levels)))
    for (j in seq_along(cols)) ipd$data[[cols[j]]] <- as.numeric(ipd$data$S_E == levels[j + 1L])
    model <- "saturated"; label <- "saturated EM total score"
  }
  list(data = ipd, columns = cols, levels = levels, model = model, label = label)
}
rd_marginal_bc_from_fit <- function(fit, y, x, wfit, wtarget, run_exact,
    seed = NULL, bc_control = list()) {
  empty <- list(estimate = NA_real_, low = NA_real_, up = NA_real_, p = NA_real_,
    status = "fit_failed", interval.available = FALSE, p.status = "fit_failed")
  if (is.null(fit) || !isTRUE(fit$converged)) return(empty)
  target <- marginal_rd_target(fit$point.est, fit$va, fit$vb, wtarget, TRUE)
  empty$estimate <- target$value
  if (!run_exact) {
    empty$status <- empty$p.status <- "BC_not_run"
    return(empty)
  }
  args <- list(param = "RD", y = y, x = x, va = fit$va, vb = fit$vb, weight = wfit,
    pars = fit$point.est, standardization.weight = wtarget, seed = seed)
  variance <- drop(crossprod(target$gradient, fit$cov %*% target$gradient))
  if (is.finite(variance) && variance > 0) args$vcov.pars <- fit$cov else args$psi.se <- .1
  tryCatch(do.call(exact_marginal_rd_target, c(args, bc_control)), error = function(e) {
    empty$status <- empty$p.status <- if (inherits(e, "marginal_exact_timeout")) "timeout" else "error"
    empty$message <- conditionMessage(e); empty
  })
}

get.compare.anchored.rd <- function(data.ipd, data.ad, param = "RD", run_exact = FALSE,
    seed = NULL, exact_seed_offset = 1000000L, bc_control = list(),
    weight_covariates = paste0("z", 1:6), score_spec = NULL,
    model_covariates = NULL, effect_modifiers = NULL,
    score_standardization = c("bc_estimated", "fixed"), fit_control = list(),
    model_weighting = c("unit", "maic"), return_full = TRUE,
    estimation_mode = c("scores", "none"),
    brm_effect_model = c("linear_atanh", "saturated")) {
  if (!identical(param, "RD")) stop("This entry point supports RD only.")
  estimation_mode <- match.arg(estimation_mode)
  brm_effect_model <- match.arg(brm_effect_model)
  if (estimation_mode == "none") {
    score_spec <- make_brm_score_spec(character(0), character(0))
    model_covariates <- effect_modifiers <- character(0)
  }
  model_weighting <- match.arg(model_weighting); score_standardization <- match.arg(score_standardization)
  fit_control <- anchored_fit_control(fit_control); validate_rd_bc_control(bc_control)
  if (!is.logical(run_exact) || length(run_exact) != 1 || is.na(run_exact)) stop("run_exact must be TRUE/FALSE.")
  if (!is.null(seed)) { rd_validate_seeds(seed); rd_validate_seeds(seed + exact_seed_offset) }
  if (!is.list(data.ipd) || !is.data.frame(data.ipd$data)) stop("data.ipd must contain individual data.")
  if (is.null(score_spec)) score_spec <- make_brm_score_spec(paste0("z", 1:3), paste0("z", 4:6))
  if (score_standardization == "bc_estimated") score_spec <- brm_scores_from_bc(score_spec, data.ad)
  prepared <- prepare_brm_score_model(data.ipd, score_spec, effect_modifiers, NULL, model_covariates)
  ipd <- prepared$data; dat <- ipd$data; n <- nrow(dat); y <- dat$y; x <- dat$x
  if (n < 2 || any(!is.finite(c(y, x))) || any(!y %in% 0:1) || any(!x %in% 0:1) || length(unique(x)) != 2)
    stop("Individual data need finite binary y/x and both treatment arms.")
  if (is.null(model_covariates)) model_covariates <- prepared$model_covariates
  if (is.null(effect_modifiers)) effect_modifiers <- prepared$effect_modifiers
  vars <- c(weight_covariates, model_covariates, effect_modifiers)
  if (!is.character(vars) || anyNA(vars) || any(!vars %in% names(dat)) || any(vars %in% c("y", "x")))
    stop("Analysis covariates must be existing baseline columns.")
  for (v in list(weight_covariates, model_covariates, effect_modifiers)) if (anyDuplicated(v)) stop("Duplicate covariates.")
  if (any(!effect_modifiers %in% model_covariates)) stop("effect_modifiers must be included in model_covariates.")
  data.ad <- brm_score_summary(data.ad, score_spec, intersect(weight_covariates, c("h_E", "h_P", "S_E", "S_P")))
  Z <- as.matrix(dat[, weight_covariates, drop = FALSE])
  mu <- select_weight_summary(data.ad$anchored.mean.x, weight_covariates, covariate_names(dat))
  Vmu <- select_weight_summary(data.ad$anchored.cov.x, weight_covariates, covariate_names(dat)) / data.ad$anchored.n
  wf <- tryCatch(get.weight(Z, mu), maic_calibration_error = function(e) e$result)
  wfit <- if (model_weighting == "unit") rep(1, n) else wf$weight
  rows <- rd_empty_rows(model_weighting); fits <- list(); bc <- list(); warnings <- list(); errors <- list()
  basis <- rd_brm_basis(ipd, score_spec, brm_effect_model)
  rows$brm_effect_model <- basis$model
  variables <- list(effect_scale = "RD", comparison_type = "anchored", method_naming_version = 1L,
    weight_covariates = weight_covariates, model_covariates = model_covariates,
    effect_modifiers = effect_modifiers, brm_effect_modifiers = basis$columns,
    brm_effect_model = basis$model, brm_effect_model_requested = brm_effect_model,
    brm_effect_basis = basis$label,
    brm_score_levels = basis$levels,
    probability_mapping = .anchored_rd_solver_env$probability_mapping,
    estimation_mode = estimation_mode,
    score_spec = score_spec, score_standardization = score_standardization,
    model_weighting = model_weighting, standardization_weighting = "maic",
    fit_control = fit_control, bc_control = bc_control, run_exact = run_exact)
  assemble <- function() structure(list(estimates = rows, marginal_estimates = rows,
    coefficient_estimates = rd_coefficient_rows(fits), fits = fits, bc = bc, weights = wf,
    model_weights = list(weight = wfit, status = if (model_weighting == "unit") "available" else wf$status),
    variables = variables, warnings = warnings, errors = errors,
    data = list(data.IPD = ipd, data.AD = data.ad)), class = c("anchored_rd_comparison", "list"))
  if (!identical(wf$status, "available")) {
    rows$status <- rows$se_status <- rows$interval_status <- rows$p_status <- "weighting_failed"
    rows$message <- wf$diagnostics$reason
    return(if (return_full) assemble() else rows)
  }
  capture <- function(label, expression) {
    msgs <- character(0)
    ans <- tryCatch(withCallingHandlers(expression, warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) { errors[[label]] <<- conditionMessage(e); NULL })
    warnings[[label]] <<- msgs; ans
  }
  ## Always use the shared fixed-weight fit covariance branch here. The full
  ## joint marginal variance (including target standardization) is computed below.
  for (label in c("brm", "brm-FC")) fits[[label]] <- capture(label,
    rd_get_estimate("AD", basis$data, wfit, gamma.ipd = numeric(0), target_mean = numeric(0),
      method = if (label == "brm") "brm" else "firth", param = "RD",
      effect_modifiers = basis$columns, non_effect_modifiers = character(0),
      model_covariates = basis$columns, weight_covariates = character(0), op_model = "marginal",
      fit_control = fit_control))
  gd <- rd_glm_design(dat, model_covariates, effect_modifiers)
  fits$IB <- capture("IB", rd_fit_identity(dat, wfit, gd, fit_control$glm_maxit))
  fits$RLPM <- capture("RLPM", rd_fit_lpm(dat, wfit, gd))
  arms <- capture("MN_eff", rd_weighted_arms(y, x, wf$weight))
  fits$MN_eff <- if (is.null(arms)) NULL else capture("MN_eff",
    rd_mn_interval(arms$p1, arms$n1, arms$p0, arms$n0))
  fits$RD_MAIC_sandwich <- capture("RD_MAIC_sandwich", rd_weighted_sandwich(y, x, wf$weight,
    Z, mu, Vmu, joint = model_weighting == "maic"))
  for (i in seq_len(nrow(rows))) {
    label <- rows$method[i]; method <- rows$fit_method[i]; fit <- fits[[method]]
    if (is.null(fit) || (method %in% c("brm", "brm-FC") && !isTRUE(fit$converged))) {
      st <- if (is.null(fit)) "fit_failed" else "fit_not_converged"
      rows$status[i] <- rows$se_status[i] <- rows$interval_status[i] <- rows$p_status[i] <- st
      rows$message[i] <- if (is.null(errors[[method]])) paste(warnings[[method]], collapse = " | ") else errors[[method]]
      next
    }
    if (rows$ci_type[i] == "BC") {
      cseed <- if (is.null(seed)) NULL else seed + exact_seed_offset
      z <- rd_marginal_bc_from_fit(fit, y, x, wfit, wf$weight, run_exact, cseed, bc_control)
      bc[[if (method == "brm") "brm.bc" else "firth.bc"]] <- z
      rows$estimate[i] <- z$estimate; rows$lower[i] <- z$low; rows$upper[i] <- z$up; rows$p.value[i] <- z$p
      rows$status[i] <- if (is.finite(z$estimate)) "available" else "target_failed"
      rows$se_status[i] <- "not_provided"
      rows$interval_status[i] <- if (z$status == "boundary_truncated") "BC_boundary_truncated" else z$status
      rows$p_status[i] <- z$p.status
      if (!is.null(z$message)) rows$message[i] <- z$message
      next
    }
    z <- capture(label, {
      if (method %in% c("MN_eff", "RD_MAIC_sandwich")) fit else {
        terms <- if (method %in% c("brm", "brm-FC")) rd_brm_terms(fit, y, x, wfit) else
          rd_regression_terms(fit, robust = label %in% c("RIB", "RLPM"))
        rd_delta(terms, wf$weight, wfit, Z, mu, Vmu,
          joint = model_weighting == "maic" && label != "IB")
      }
    })
    if (is.null(z)) {
      rows$status[i] <- rows$se_status[i] <- rows$interval_status[i] <- rows$p_status[i] <- "standardization_failed"
      rows$message[i] <- errors[[label]]
      ## Recover an independently computable point after a variance failure.
      point <- tryCatch(if (method %in% c("brm", "brm-FC"))
        marginal_rd_target(fit$point.est, fit$va, fit$vb, wf$weight) else
        weighted.mean(drop((fit$design$X1 - fit$design$X0) %*% fit$coefficients), wf$weight), error = function(e) NA_real_)
      if (is.finite(point)) { rows$estimate[i] <- point; rows$status[i] <- "available" }
      next
    }
    for (key in intersect(names(z), names(rows))) rows[i, key] <- z[[key]]
    rows$status[i] <- if (is.finite(z$estimate)) "available" else "target_failed"
    rows$se_status[i] <- if (is.finite(z$se) && z$se >= 0) "available" else "not_provided"
    rows$interval_status[i] <- if (all(is.finite(c(z$lower, z$upper)))) "available" else "se_failed"
    rows$p_status[i] <- if (is.finite(z$p.value)) "available" else "p_failed"
    if (rows$method[i] == "MN_eff") rows$se_status[i] <- "diagnostic_ci_width_only"
    if (label == "IB" && model_weighting == "maic") rows$uncertainty[i] <- "fixed MAIC fitting/standardization weights; model-based"
    rows$message[i] <- paste(unique(warnings[[label]]), collapse = " | ")
  }
  out <- assemble(); if (return_full) out else rows
}

rd_coefficient_rows <- function(fits) {
  rows <- list()
  for (method in names(fits)) {
    fit <- fits[[method]]
    if (method %in% c("brm", "brm-FC") && !is.null(fit$point.est)) {
      pa <- ncol(fit$va); cf <- fit$coefficients
      rows[[method]] <- data.frame(method = method, component = c(rep("alpha", pa), rep("beta", ncol(fit$vb))),
        term = c(colnames(fit$va), colnames(fit$vb)), parameter = c(rep("atanh_RD", pa), rep("log_OP", ncol(fit$vb))),
        estimate = fit$point.est, se = fit$se.est, lower = fit$conf.lower, upper = fit$conf.upper,
        status = if (isTRUE(fit$converged)) "available" else "fit_not_converged", row.names = NULL)
    } else if (method %in% c("IB", "RLPM") && !is.null(fit$coefficients))
      rows[[method]] <- data.frame(method = method, component = "regression", term = names(fit$coefficients),
        parameter = "identity_coefficient", estimate = fit$coefficients, se = sqrt(diag(fit$covariance)),
        lower = NA_real_, upper = NA_real_, status = "available", row.names = NULL)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}
run.anchored.rd <- function(design, n_AC, n_BC = n_AC, seed = NULL,
    return_data = TRUE, return_full = TRUE, ...) {
  rd_validate_sizes(n_AC, n_BC)
  if (!is.list(design) || !identical(design$effect_scale, "RD")) stop("An RD design is required.")
  if (!is.null(seed)) { rd_validate_seeds(seed); set.seed(seed) }
  sim <- generate_anchored_data(n_AC, n_BC, design)
  result <- get.compare.anchored.rd(sim$data.IPD, sim$data.AD, seed = seed, return_full = TRUE, ...)
  result$simulation <- list(design = design, n_AC = n_AC, n_BC = n_BC, seed = seed)
  result$estimates$population_RD_AC_in_BC <- unname(design$marginal_truth["BC", "RD_AC"])
  result$estimates$population_RD_BC_in_BC <- unname(design$marginal_truth["BC", "RD_BC"])
  result$estimates$population_RD_AB_in_BC <- unname(design$marginal_truth["BC", "RD_AB"])
  result$marginal_estimates <- result$estimates
  if (return_data) { result$data$data.AD.or <- sim$data.AD.or; result$data$design <- design }
  if (return_full) result else result$estimates
}
