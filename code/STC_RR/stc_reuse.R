# STC post-processing of saved RR v2 score_unit fits. Base R only.
# No generators, optimizers, model fitting, or simulation inversion are loaded.
stc_version <- "STC_RR_reuse_20261026_v1"
stc_methods <- c("brm", "brm-FC", "LB", "LP", "RLB", "RLP",
                 "LB-FC", "RLB-FC", "LP-FC", "RLP-FC", "LP-BFC", "RLP-BFC")
stc_scopes <- c("fixed_bc_summary", "bc_mean_delta")
stc_cols <- paste0("z", 1:6)

stc_get <- function(x, name) if (is.list(x)) x[[name, exact = TRUE]] else NULL
stc_stop <- function(status, message) {
  stop(structure(list(message = message, call = NULL, status = status),
                 class = c("stc_error", "error", "condition")))
}
stc_require <- function(ok, status, message) {
  if (!isTRUE(ok)) stc_stop(status, message)
}
stc_error_status <- function(e, default = "postprocessing_failed") {
  if (inherits(e, "stc_error")) e$status else default
}

stc_scenarios <- function() {
  x <- data.frame(scenario_id = 1:8, event = rep(c("common", "rare"), each = 4),
    n = rep(c(100L, 100L, 500L, 500L), 2), ess_ratio = rep(c(.15, .5), 4))
  x$cell <- sprintf("S%02d_%s_n%d_ess%03d_score_unit", x$scenario_id,
                    x$event, x$n, as.integer(round(100 * x$ess_ratio)))
  x
}

# The reported BC covariance is used ONLY for mean-delta uncertainty here.
# Independence is an explicit target reconstruction assumption, not inferred
# from a zero sample covariance and not supplied by the DGP's true probabilities.
stc_binary_target <- function(ad) {
  mu <- stc_get(ad, "anchored.mean.x")
  stc_require(is.numeric(mu) && !is.null(names(mu)) && !anyDuplicated(names(mu)) &&
    all(stc_cols %in% names(mu)), "bc_mean_missing", "Six named BC means are required.")
  mu <- mu[stc_cols]
  stc_require(all(is.finite(mu)) && all(mu >= 0 & mu <= 1),
              "bc_mean_invalid", "BC means must be finite probabilities in [0,1].")
  z <- expand.grid(setNames(rep(list(0:1), 6), stc_cols))
  factors <- vapply(1:6, function(k) ifelse(z[[k]] == 1, mu[k], 1 - mu[k]), numeric(64))
  omega <- apply(factors, 1, prod)
  D <- vapply(1:6, function(k) (2 * z[[k]] - 1) *
    apply(factors[, -k, drop = FALSE], 1, prod), numeric(64))
  colnames(D) <- stc_cols
  stc_require(abs(sum(omega) - 1) < 1e-12, "target_mass_invalid", "Target mass does not sum to one.")
  list(data = z, mu = mu, weight = omega, Dweight = D,
       boundary = any(mu == 0 | mu == 1), ad = ad)
}

# Reapply the SAVED affine score coding; do not recompute it from true moments.
stc_apply_scores <- function(data, spec) {
  stc_require(inherits(spec, "brm_score_spec"), "score_spec_missing", "Saved brm_score_spec is required.")
  for (group in c("E", "P")) {
    key <- paste0("h_", group)
    cols <- spec[[if (group == "E") "em_covariates" else "non_em_covariates"]]
    stc_require(is.character(cols) && all(cols %in% names(data)),
                "score_spec_invalid", "Score components are missing from prediction data.")
    coeff <- spec$coefficients[cols]
    stc_require(length(coeff) == length(cols) && all(is.finite(coeff)) &&
      length(spec$scale[key]) == 1 && is.finite(spec$scale[key]) && spec$scale[key] > 0 &&
      length(spec$center[key]) == 1 && is.finite(spec$center[key]),
      "score_spec_invalid", "Saved score coefficients, center or scale are invalid.")
    total <- numeric(nrow(data))
    for (k in cols) {
      v <- data[[k]]
      stc_require(is.numeric(v) && all(is.finite(v)) && all(v %in% 0:1),
                  "score_data_invalid", "Raw score components must be binary.")
      if (k %in% spec$reverse_covariates) v <- 1 - v
      total <- total + spec$coefficients[k] * v
    }
    data[[paste0("S_", group)]] <- total
    data[[key]] <- (total - spec$center[key]) / spec$scale[key]
  }
  data
}

stc_validate_record <- function(record, expected = NULL) {
  s <- stc_get(record, "settings"); old <- stc_get(record, "summary")
  stc_require(is.list(s), "settings_missing", "Saved settings are missing.")
  stc_require(identical(s$dgm, "webster2024_rr_ess_em_v2") &&
    identical(as.integer(s$dgm_scenario$version), 2L) &&
    identical(s$dgm_scenario$risk_target_version, "BC_CAB_20261015"),
    "dgm_mismatch", "Expected RR v2 / BC_CAB_20261015.")
  stc_require(identical(s$model_weighting, "unit") && identical(s$ess_source, "diffuse") &&
    identical(s$em_level, "high") && identical(s$weight_mode, "all"),
    "analysis_mismatch", "Only diffuse/high/all score_unit fits may be reused.")
  mode <- s$estimation_mode
  if (is.null(mode) && is.data.frame(old) && "estimation_mode" %in% names(old))
    mode <- unique(old$estimation_mode)
  stc_require(identical(mode, "scores"), "analysis_mismatch", "Saved estimation_mode must be scores.")
  grid <- stc_scenarios()
  i <- which(grid$event == s$event & grid$n == s$n_AC & grid$n == s$n_BC &
               abs(grid$ess_ratio - s$ess_ratio) < 1e-12)
  stc_require(length(i) == 1L, "scenario_mismatch", "Record is not one of the eight requested scenarios.")
  stc_require(length(s$seed) == 1L && is.finite(s$seed) && s$seed %in% 1:1000,
              "seed_invalid", "Expected an original data seed in 1..1000.")
  if (!is.null(expected)) {
    stc_require(identical(as.integer(grid$scenario_id[i]), as.integer(expected$scenario_id)) &&
      identical(as.integer(s$seed), as.integer(expected$seed)),
      "source_identity_mismatch", "RDS settings disagree with source cell/seed.")
  }
  if (is.data.frame(old) && nrow(old)) {
    for (key in intersect(c("event", "n_AC", "n_BC", "ess_ratio", "seed", "dgm", "model_weighting"), names(old)))
      stc_require(all(!is.na(old[[key]]) & old[[key]] == s[[key]]),
                  "source_identity_mismatch", paste("Summary/settings disagree:", key))
  }
  list(settings = s, scenario = grid[i, , drop = FALSE])
}

# Stable RR/OP conversion, frozen from the original RR helper (2026-10-24).
# pmin(logp,0) below is the original roundoff guard in the BRM algebra;
# it is not probability truncation applied to GLM predictions.
stc_logadd <- function(a, b) {
  m <- pmax(a, b); ans <- m + log1p(exp(-abs(a - b)))
  ans[is.infinite(m) & m < 0] <- -Inf
  ans
}
stc_brm_risks <- function(logrr, logop) {
  lr1 <- ifelse(logrr > 0, logrr + log(-expm1(-pmax(logrr, 0))),
                log(-expm1(pmin(logrr, 0))))
  lb <- logop + stc_logadd(0, logrr)
  ld <- .5 * stc_logadd(2 * logop + 2 * lr1, log(4) + logrr + logop)
  l0 <- log(2) + logop - stc_logadd(lb, ld); l1 <- logrr + l0
  l0 <- pmin(l0, 0); l1 <- pmin(l1, 0)
  q0 <- -expm1(l0); q1 <- -expm1(l1); den <- q0 + q1
  list(p0 = exp(l0), p1 = exp(l1), t0 = -q0 / den,
       t1 = q1 / den, s = q0 * q1 / den)
}

stc_unwrap <- function(fit) {
  inner <- stc_get(fit, "fit") # exact lookup: glm$fit would match fitted.values!
  if (is.null(inner)) fit else inner
}
stc_check_fit <- function(fit, brm = FALSE) {
  stc_require(is.list(fit), "fit_missing", "Saved model fit is missing.")
  actual <- if (brm) fit else stc_unwrap(fit)
  stc_require(is.list(actual), "fit_schema_invalid", "Saved model object is invalid.")
  stc_require(isTRUE(actual$converged) && !isFALSE(stc_get(fit, "converged")),
              "fit_not_converged", "Original model did not converge; no refit attempted.")
  actual
}
stc_covariance <- function(fit, method) {
  if (method %in% c("brm", "brm-FC")) return(stc_get(fit, "cov"))
  # A missing robust/FC wrapper covariance must NEVER silently fall back to
  # a conventional model-based covariance from the underlying glm.
  if (!method %in% c("LB", "LP")) return(stc_get(fit, "vcov"))
  actual <- stc_unwrap(fit)
  # Reconstruct a standard glm covariance from the saved QR/dispersion only.
  # stats::vcov does not refit a model.
  if (inherits(actual, "glm")) return(stats::vcov(actual))
  NULL
}
stc_glm_design <- function(actual, nd) {
  beta <- stc_get(actual, "coefficients")
  stc_require(is.numeric(beta) && length(beta) > 0 && !is.null(names(beta)) &&
    !anyDuplicated(names(beta)) && all(is.finite(beta)),
    "coefficients_invalid", "GLM coefficients are missing, aliased or nonfinite.")
  tt <- stc_get(actual, "terms")
  stc_require(inherits(tt, "terms"), "model_terms_missing", "Saved GLM terms are required.")
  tt <- stats::delete.response(tt)
  matrices <- lapply(0:1, function(arm) {
    d <- nd; d$x <- arm
    frame <- stats::model.frame(tt, d, xlev = stc_get(actual, "xlevels"), na.action = stats::na.fail)
    X <- stats::model.matrix(tt, frame, contrasts.arg = stc_get(actual, "contrasts"))
    stc_require(setequal(colnames(X), names(beta)), "model_terms_mismatch", "GLM prediction columns disagree with coefficients.")
    X[, names(beta), drop = FALSE]
  })
  list(beta = beta, X0 = matrices[[1]], X1 = matrices[[2]])
}
stc_glm_terms <- function(actual, nd, beta = NULL) {
  d <- stc_glm_design(actual, nd)
  if (is.null(beta)) beta <- d$beta
  stc_require(length(beta) == ncol(d$X0), "coefficients_invalid", "Wrong GLM parameter length.")
  family <- stc_get(actual, "family")
  stc_require(is.list(family) && identical(family$link, "log"),
              "model_link_mismatch", "Only the saved log-link models are supported.")
  p0 <- exp(drop(d$X0 %*% beta)); p1 <- exp(drop(d$X1 %*% beta))
  list(p0 = p0, p1 = p1, dp0 = d$X0 * p0, dp1 = d$X1 * p1,
       parameter_names = names(d$beta), parameter_count = length(beta))
}

stc_brm_design <- function(fit, record, target, spec) {
  stc_require(identical(fit$op_spec$mode, "marginal") &&
    identical(fit$op_spec$conditioning_covariates, "h_E"),
    "brm_design_unsupported", "Expected the original reduced BRM with saturated EM-score OP.")
  ac <- stc_get(stc_get(stc_get(record$result, "data"), "data.IPD"), "data")
  stc_require(is.data.frame(ac), "ac_design_missing", "Saved AC design rows are required for BRM mapping.")
  observed <- stc_apply_scores(ac[, stc_cols, drop = FALSE], spec)
  nd <- stc_apply_scores(target$data, spec)
  va <- fit$va; vb <- fit$vb
  stc_require(is.matrix(va) && is.matrix(vb) && nrow(va) == nrow(ac) &&
    nrow(vb) == nrow(ac) && ncol(va) == 2L && all(is.finite(c(va, vb))) &&
    max(abs(va[, 1] - 1)) < 1e-10 && max(abs(va[, 2] - observed$h_E)) < 1e-8,
    "brm_design_mismatch", "Saved BRM design does not match the saved score coding/AC row order.")
  # Retain the original cell coding/order, including its nuisance reference.
  map <- match(nd$S_E, observed$S_E)
  missing <- is.na(map) & target$weight > 0
  stc_require(!any(missing), "target_support_missing", paste0(
    "No saved nuisance cell for target EM score(s): ", paste(unique(nd$S_E[missing]), collapse = ",")))
  for (score in unique(observed$S_E)) {
    ii <- which(observed$S_E == score)
    stc_require(max(abs(sweep(vb[ii, , drop = FALSE], 2, vb[ii[1], ], "-"))) < 1e-10,
      "brm_design_mismatch", "Nuisance design is not constant within an EM score.")
  }
  # Unobserved ZERO-mass cells cannot affect the fixed target. Mean-delta
  # inference is separately refused for boundary target means.
  map[is.na(map)] <- 1L
  list(va = va[map, , drop = FALSE], vb = vb[map, , drop = FALSE])
}
stc_brm_terms <- function(fit, design, beta = fit$point.est) {
  pa <- ncol(design$va); pb <- ncol(design$vb)
  stc_require(is.numeric(beta) && length(beta) == pa + pb && all(is.finite(beta)),
              "coefficients_invalid", "BRM full coefficient vector is invalid.")
  z <- stc_brm_risks(drop(design$va %*% beta[seq_len(pa)]),
                    drop(design$vb %*% beta[pa + seq_len(pb)]))
  list(p0 = z$p0, p1 = z$p1,
    dp0 = cbind(design$va * (z$p0 * z$t0), design$vb * (z$p0 * z$s)),
    dp1 = cbind(design$va * (z$p1 * z$t1), design$vb * (z$p1 * z$s)),
    parameter_names = names(beta), parameter_count = length(beta))
}

stc_evaluate <- function(terms, target) {
  w <- target$weight; active <- w > 0
  stc_require(length(terms$p0) == length(w) && length(terms$p1) == length(w) &&
    all(is.finite(c(terms$p0[active], terms$p1[active]))),
    "prediction_invalid", "Target predictions are nonfinite.")
  m0 <- sum(w[active] * terms$p0[active]); m1 <- sum(w[active] * terms$p1[active])
  stc_require(all(is.finite(c(m0, m1))) && min(m0, m1) > 0,
              "standardized_risk_invalid", "Standardized risks must be finite and positive.")
  # Do not require each GLM prediction <=1 and do not truncate any prediction.
  gradient <- colSums(w[active] * (terms$dp1[active, , drop = FALSE] / m1 -
                                  terms$dp0[active, , drop = FALSE] / m0))
  gmu <- colSums(target$Dweight * (terms$p1 / m1 - terms$p0 / m0))
  pred <- c(terms$p0[active], terms$p1[active])
  list(estimate = log(m1) - log(m0), risk_C = m0, risk_A = m1,
    gradient = gradient, gradient_mu = gmu,
    prediction_out_of_bounds = any(pred < 0 | pred > 1),
    n_predictions_out_of_bounds = sum(pred < 0 | pred > 1),
    min_prediction = min(pred), max_prediction = max(pred))
}
stc_variance <- function(g, V, parameter_names = NULL) {
  stc_require(is.numeric(g) && all(is.finite(g)) && is.matrix(V) &&
    identical(dim(V), c(length(g), length(g))) && all(is.finite(V)),
    "covariance_invalid", "Complete finite parameter covariance matrix is required.")
  if (!is.null(parameter_names) && !is.null(rownames(V)) && !is.null(colnames(V))) {
    stc_require(setequal(rownames(V), parameter_names) && setequal(colnames(V), parameter_names),
                "covariance_invalid", "Covariance parameter names do not match coefficients.")
    V <- V[parameter_names, parameter_names, drop = FALSE]
  }
  scale <- max(1, max(abs(V)))
  stc_require(max(abs(V - t(V))) <= 1e-8 * scale,
              "covariance_invalid", "Covariance matrix is not symmetric.")
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  stc_require(min(ev) >= -1e-8 * scale, "covariance_invalid", "Covariance matrix is not positive semidefinite.")
  v <- drop(crossprod(g, V %*% g))
  stc_require(is.finite(v) && v >= -1e-10 * scale * max(1, sum(g^2)),
              "covariance_invalid", "Delta variance is negative or nonfinite.")
  max(v, 0) # roundoff only; no modification to saved covariance
}
stc_target_variance <- function(point, target) {
  stc_require(!target$boundary, "target_mean_boundary", "BC-mean delta interval is not regular at a boundary mean.")
  ad <- target$ad; V <- stc_get(ad, "anchored.cov.x"); n <- stc_get(ad, "anchored.n")
  stc_require(is.numeric(n) && length(n) == 1 && is.finite(n) && n > 1 &&
    is.matrix(V) && all(stc_cols %in% rownames(V)) && all(stc_cols %in% colnames(V)),
    "bc_covariance_missing", "Reported BC covariance and sample size are required for mean-delta inference.")
  stc_variance(point$gradient_mu, V[stc_cols, stc_cols, drop = FALSE] / n, stc_cols)
}

stc_empty_rows <- function(scenario, seed, status = "pending", message = "", level = .95) {
  rows <- expand.grid(uncertainty = stc_scopes, method = stc_methods, stringsAsFactors = FALSE)
  rows$scenario_id <- scenario$scenario_id; rows$cell <- scenario$cell
  rows$event <- scenario$event; rows$n_AC <- rows$n_BC <- scenario$n
  rows$ess_ratio <- scenario$ess_ratio; rows$seed <- as.integer(seed)
  rows$target <- "A/C marginal log RR in BC target"
  rows$dgm <- "webster2024_rr_ess_em_v2"; rows$analysis <- "STC_score_unit_reuse"
  rows$model_weighting <- "unit"; rows$standardization <- "independent_binary_BC"
  rows$score_coding <- "saved_fixed"
  rows$risk_interpretation <- ifelse(rows$method %in% c("brm", "brm-FC"),
    "AC reduced-model risks averaged over BC EM distribution; not BC absolute risks",
    "model risks averaged over reconstructed BC covariate distribution")
  rows$ci_type <- "Wald"; rows$level <- level; rows$version <- stc_version
  for (key in c("estimate", "se", "low", "up", "p", "RR", "RR_low", "RR_up", "risk_C", "risk_A",
                "variance_fit", "variance_target", "population_log_RR_AC_in_BC", "min_prediction", "max_prediction"))
    rows[[key]] <- NA_real_
  rows$prediction_out_of_bounds <- NA; rows$n_predictions_out_of_bounds <- NA_integer_
  rows$status <- rows$point_status <- rows$se_status <- status
  rows$message <- message; rows$source_file <- ""; rows$source_md5 <- ""
  rows
}
stc_reference_truth <- function(record) {
  x <- stc_get(record, "summary")
  if (!is.data.frame(x) || !"population_log_RR_AC_in_BC" %in% names(x)) return(NA_real_)
  v <- unique(x$population_log_RR_AC_in_BC[is.finite(x$population_log_RR_AC_in_BC)])
  if (length(v) != 1L) return(NA_real_)
  v
}
stc_pair_rows <- function(rows, record = NULL) {
  out <- rows
  old <- stc_get(record, "summary")
  for (key in c("estimate", "se", "low", "up")) out[[paste0("maic_", key)]] <- NA_real_
  out$maic_status <- "source_unavailable"; out$maic_uncertainty <- ""
  if (is.data.frame(old) && all(c("method", "ci_type", "estimate", "se", "low", "up", "status") %in% names(old))) {
    old <- old[old$ci_type == "Wald", , drop = FALSE]
    for (m in stc_methods) {
      ii <- which(out$method == m); z <- old[old$method == m, , drop = FALSE]
      if (nrow(z) != 1L) next
      for (key in c("estimate", "se", "low", "up")) out[ii, paste0("maic_", key)] <- z[[key]][1]
      out$maic_status[ii] <- z$status[1]
      if ("uncertainty" %in% names(z)) out$maic_uncertainty[ii] <- z$uncertainty[1]
    }
  }
  out$stc_point_valid <- is.finite(out$estimate) & out$point_status == "available"
  out$maic_point_valid <- is.finite(out$maic_estimate) & out$maic_status %in% c("available", "se_failed")
  out$paired_point_valid <- out$stc_point_valid & out$maic_point_valid
  out$stc_interval_valid <- out$status == "available" & is.finite(out$low) & is.finite(out$up)
  out$maic_interval_valid <- out$maic_status == "available" & is.finite(out$maic_low) & is.finite(out$maic_up)
  out$paired_interval_valid <- out$stc_interval_valid & out$maic_interval_valid
  out$stc_minus_maic <- ifelse(out$paired_point_valid, out$estimate - out$maic_estimate, NA_real_)
  out
}

# Public interface. Returns summary (24 rows), paired (24 rows), diagnostics.
# No field in record is modified. Raw BC individual data and oracle risks are
# never accessed; reference truth is read only for subsequent error evaluation.
stc_from_record <- function(record, expected = NULL, level = .95, source_file = "") {
  stc_require(length(level) == 1 && is.finite(level) && level > 0 && level < 1,
              "level_invalid", "level must be in (0,1).")
  checked <- stc_validate_record(record, expected)
  rows <- stc_empty_rows(checked$scenario, checked$settings$seed, level = level)
  rows$source_file <- source_file
  rows$population_log_RR_AC_in_BC <- stc_reference_truth(record) # evaluation only
  if (!isTRUE(record$complete)) {
    rows$status <- rows$point_status <- rows$se_status <- "source_incomplete"
    rows$message <- "Original seed record is incomplete; no fitting or simulation attempted."
    return(list(summary = rows, paired = stc_pair_rows(rows, record), diagnostics = data.frame()))
  }
  result <- stc_get(record, "result")
  stc_require(is.list(result), "result_missing", "Saved result object is missing.")
  stc_require(is.list(stc_get(result, "fits")), "fits_missing", "Saved fit collection is missing.")
  fit_weights <- stc_get(stc_get(result, "model_weights"), "weight")
  stc_require(identical(stc_get(stc_get(result, "variables"), "model_weighting"), "unit") &&
    is.numeric(fit_weights) && length(fit_weights) == checked$settings$n_AC &&
    all(is.finite(fit_weights)) && max(abs(fit_weights - 1)) < 1e-12,
    "analysis_mismatch", "Saved actual model fitting weights must be all ones.")
  spec <- stc_get(stc_get(result, "variables"), "score_spec")
  stc_require(inherits(spec, "brm_score_spec") &&
    identical(spec$em_covariates, stc_cols[1:3]) && identical(spec$non_em_covariates, stc_cols[4:6]),
    "score_spec_missing", "Original six-covariate score specification is required.")
  ad <- stc_get(stc_get(result, "data"), "data.AD")
  if (!is.null(stc_get(ad, "anchored.n")))
    stc_require(identical(as.integer(ad$anchored.n), as.integer(checked$settings$n_BC)),
      "source_identity_mismatch", "Reported BC sample size disagrees with saved settings.")
  target <- stc_binary_target(ad)
  nd <- stc_apply_scores(target$data, spec)
  diagnostics <- list()
  for (method in stc_methods) {
    ii <- which(rows$method == method); fit <- stc_get(stc_get(result, "fits"), method)
    point <- tryCatch({
      brm <- method %in% c("brm", "brm-FC")
      actual <- stc_check_fit(fit, brm)
      terms <- if (brm) stc_brm_terms(fit, stc_brm_design(fit, record, target, spec)) else
        stc_glm_terms(actual, nd)
      stc_evaluate(terms, target)
    }, error = identity)
    if (inherits(point, "error")) {
      rows$status[ii] <- rows$point_status[ii] <- rows$se_status[ii] <- stc_error_status(point)
      rows$message[ii] <- conditionMessage(point)
      next
    }
    for (key in intersect(names(point), names(rows))) rows[ii, key] <- point[[key]]
    rows$point_status[ii] <- "available"; rows$RR[ii] <- exp(point$estimate)
    vfit <- tryCatch(stc_variance(point$gradient, stc_covariance(fit, method), terms$parameter_names), error = identity)
    vtgt <- tryCatch(stc_target_variance(point, target), error = identity)
    for (j in ii) {
      if (inherits(vfit, "error")) {
        rows$status[j] <- rows$se_status[j] <- "se_failed"
        rows$message[j] <- paste(stc_error_status(vfit), conditionMessage(vfit), sep = ": ")
        next
      }
      rows$variance_fit[j] <- vfit
      if (rows$uncertainty[j] == "bc_mean_delta" && inherits(vtgt, "error")) {
        rows$status[j] <- rows$se_status[j] <- stc_error_status(vtgt)
        rows$message[j] <- conditionMessage(vtgt)
        next
      }
      extra <- if (rows$uncertainty[j] == "bc_mean_delta") vtgt else 0
      rows$variance_target[j] <- extra; variance <- vfit + extra
      if (!is.finite(variance) || variance <= 0) {
        rows$status[j] <- rows$se_status[j] <- "se_failed"
        rows$message[j] <- "Nonpositive or nonfinite total variance."
        next
      }
      se <- sqrt(variance); q <- stats::qnorm((1 + level) / 2)
      rows$se[j] <- se; rows$low[j] <- point$estimate - q * se
      rows$up[j] <- point$estimate + q * se
      rows$p[j] <- 2 * stats::pnorm(-abs(point$estimate / se))
      rows$RR_low[j] <- exp(rows$low[j]); rows$RR_up[j] <- exp(rows$up[j])
      rows$status[j] <- rows$se_status[j] <- "available"
      if (point$prediction_out_of_bounds) rows$message[j] <- "Untruncated GLM prediction(s) outside [0,1]; point retained."
    }
    diagnostics[[method]] <- data.frame(method = method, target_mass = sum(target$weight),
      parameter_count = length(point$gradient), target_mean_boundary = target$boundary,
      prediction_out_of_bounds = point$prediction_out_of_bounds,
      reduced_brm = method %in% c("brm", "brm-FC"), stringsAsFactors = FALSE)
  }
  list(summary = rows, paired = stc_pair_rows(rows, record),
       diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame())
}

stc_finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
stc_metrics <- function(estimate, low, up, valid_point, valid_interval, truth, planned) {
  point <- valid_point & is.finite(estimate); interval <- valid_interval & is.finite(low) & is.finite(up)
  error <- estimate[point] - truth
  covered <- interval & is.finite(truth) & low <= truth & up >= truth
  n <- sum(interval); coverage <- if (n && is.finite(truth)) sum(covered) / n else NA_real_
  data.frame(planned = planned, points_returned = sum(point), intervals_returned = n,
    point_return_rate = sum(point) / planned, interval_return_rate = n / planned,
    truth = truth, bias = stc_finite_mean(error),
    empirical_sd = if (sum(point) > 1L) stats::sd(estimate[point]) else NA_real_,
    rmse = if (any(is.finite(error))) sqrt(mean(error[is.finite(error)]^2)) else NA_real_,
    coverage_valid = coverage,
    coverage_mcse = if (is.finite(coverage)) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    returned_and_covered_rate = if (is.finite(truth)) sum(covered) / planned else NA_real_,
    mean_width = stc_finite_mean((up - low)[interval]))
}
stc_summarize <- function(paired, planned_seeds) {
  groups <- split(seq_len(nrow(paired)), interaction(paired$scenario_id, paired$method,
                                                  paired$uncertainty, drop = TRUE))
  marginal <- common <- differences <- list()
  for (i in seq_along(groups)) {
    x <- paired[groups[[i]], , drop = FALSE]
    meta <- x[1, c("scenario_id", "event", "n_AC", "n_BC", "ess_ratio", "method", "uncertainty")]
    truth <- unique(x$population_log_RR_AC_in_BC[is.finite(x$population_log_RR_AC_in_BC)])
    truth <- if (length(truth) == 1L) truth else NA_real_
    for (path in c("STC", "MAIC_score_unit")) {
      maic <- path == "MAIC_score_unit"
      e <- if (maic) x$maic_estimate else x$estimate
      lo <- if (maic) x$maic_low else x$low; hi <- if (maic) x$maic_up else x$up
      vp <- if (maic) x$maic_point_valid else x$stc_point_valid
      vi <- if (maic) x$maic_interval_valid else x$stc_interval_valid
      k <- paste(i, path)
      marginal[[k]] <- cbind(meta, path = path, stc_metrics(e, lo, hi, vp, vi, truth, planned_seeds))
      common[[k]] <- cbind(meta, path = path, stc_metrics(e, lo, hi,
        x$paired_point_valid, x$paired_interval_valid, truth, planned_seeds))
    }
    vp <- x$paired_point_valid; vi <- x$paired_interval_valid
    dsq <- (x$estimate - truth)^2 - (x$maic_estimate - truth)^2
    differences[[i]] <- cbind(meta, planned = planned_seeds, paired_points = sum(vp),
      paired_intervals = sum(vi), mean_stc_minus_maic = stc_finite_mean(x$stc_minus_maic[vp]),
      mean_squared_error_difference = stc_finite_mean(dsq[vp]),
      mean_width_difference = stc_finite_mean(((x$up - x$low) - (x$maic_up - x$maic_low))[vi]))
  }
  list(method_summary = do.call(rbind, marginal), paired_common_summary = do.call(rbind, common),
       paired_differences = do.call(rbind, differences))
}
