## Cluster export: all methods, with marginal Wald and exact BRM intervals.
## Source after load_anchored.R/source_anchored(); this file does not launch work.
## Wald equations retain the project's archived marginal_standardization.R
## convention, including the first-order unadjusted-score convention for FC.

maic_method_plan <- function() rbind(data.frame(
  method = c("brm", "brm", "brm-FC", "brm-FC",
    "LB", "LP", "RLP", "RLB", "CMH", "CMH_sandwich", "brmad", "brmad"),
  ci_type = c(rep(c("Wald", "Exact"), 2), rep("Wald", 7), "Exact"),
  legacy_method = c("brm", "brm-BC", "brm-FC", "brm-FC-BC",
    "LB", "LP", "RLP", "RLB", "CMH", "CMH_sandwich", "brmad", "brm_b-BC"),
  stringsAsFactors = FALSE),
  data.frame(method = anchored_glm_fc_method_names(), ci_type = "Wald",
    legacy_method = anchored_glm_fc_method_names(), stringsAsFactors = FALSE))

maic_wrapped_glm <- function(method) method %in%
  c("RLP", "RLB", anchored_glm_fc_method_names())

maic_empty_rows <- function(model_weighting = c("unit", "maic")) {
  model_weighting <- match.arg(model_weighting)
  out <- maic_method_plan()
  cmh <- out$method %in% c("CMH", "CMH_sandwich")
  out$component <- "marginal"
  out$target <- "A/C marginal log RR in BC target"
  out$uncertainty <- ifelse(out$ci_type == "Exact", "fixed covariates and both sets of weights",
    "fixed fitting and standardization weights and covariates")
  out$uncertainty[cmh] <- "fixed MAIC weights; stratified count influence function"
  joint_methods <- c("brm", "brm-FC", "RLP", "RLB", "RLB-FC", "RLP-FC", "RLP-BFC")
  if (model_weighting == "maic") {
    out$uncertainty[out$ci_type == "Wald" & out$method %in% joint_methods] <-
      "joint sandwich: weights, empirical standardization and AD means"
    out$uncertainty[out$method == "CMH_sandwich"] <- "CMH joint sandwich"
  }
  out$fitting_weights <- ifelse(cmh | model_weighting == "maic", "MAIC", "ones")
  out$standardization_weights <- "MAIC"
  for (key in c("estimate", "se", "low", "up", "p", "RR", "RR_low", "RR_up",
                "risk_C", "risk_A", "variance_ipd", "variance_target", "elapsed_seconds")) out[[key]] <- NA_real_
  out$status <- out$se_status <- out$p_status <- "pending"
  out$se_status[out$ci_type == "Exact"] <- "not_provided"
  out$message <- ""
  out
}

maic_brm_terms <- function(fit, y, x, w) {
  va <- fit$va; vb <- fit$vb; pa <- ncol(va); pb <- ncol(vb); n <- length(y)
  a <- fit$point.est[seq_len(pa)]; b <- fit$point.est[pa + seq_len(pb)]
  tt <- firth_rr_derivatives(a, b, x, va, vb)
  p <- tt$p; v <- p * (1 - p)
  dp <- cbind(va * tt$d.a, vb * tt$d.b)
  score <- dp * ((y - p) / v); curv <- (y - p) / v
  hp <- rbind(cbind(crossprod(va, va * (w * curv * tt$d.aa)),
                    crossprod(va, vb * (w * curv * tt$d.ab))),
              cbind(crossprod(vb, va * (w * curv * tt$d.ab)),
                    crossprod(vb, vb * (w * curv * tt$d.bb))))
  H <- (crossprod(dp, dp * (w * (y / p^2 + (1 - y) / (1 - p)^2))) - hp) / n
  ## Use the same stable risks/gradient as exact_marginal_target for the target.
  z <- .marginal_exact_risks(drop(va %*% a), drop(vb %*% b))
  list(p0 = z$p0, p1 = z$p1,
    dp0 = cbind(va * (z$p0 * z$t0), vb * (z$p0 * z$s)),
    dp1 = cbind(va * (z$p1 * z$t1), vb * (z$p1 * z$s)),
    score = score, H = H, covariance = fit$cov)
}

maic_glm_target_terms <- function(fit, data) {
  cf <- coef(fit)
  if (any(!is.finite(cf))) stop("Nonfinite or aliased GLM coefficients.")
  d0 <- d1 <- data; d0$x <- 0; d1$x <- 1
  tt <- delete.response(terms(fit))
  X0 <- model.matrix(tt, d0, contrasts.arg = fit$contrasts)
  X1 <- model.matrix(tt, d1, contrasts.arg = fit$contrasts)
  p0 <- exp(drop(X0 %*% cf)); p1 <- exp(drop(X1 %*% cf))
  ## Export the untruncated plug-in target for every GLM. Counterfactual
  ## predictions outside [0,1] do not discard a converged, computable point.
  ## maic_target_point() still requires finite, positive standardized means;
  ## the factual fitting domain and convergence checks remain in the fitters.
  list(p0=p0,p1=p1,dp0=X0*p0,dp1=X1*p1)
}

maic_target_point <- function(terms, v) {
  v <- v / mean(v)
  m0 <- mean(v * terms$p0); m1 <- mean(v * terms$p1)
  if (any(!is.finite(c(m0,m1))) || min(m0,m1)<=0) stop("Invalid standardized risk.")
  list(estimate=log(m1)-log(m0),risk_C=m0,risk_A=m1)
}

maic_glm_terms <- function(fit, data, w) {
  n <- nrow(data); target <- maic_glm_target_terms(fit,data)
  X <- model.matrix(fit)
  p <- fitted(fit); y <- data$y
  if (inherits(fit, "glm_fc_fit")) {
    score <- fit$score_rows
    ## The FC adjustment uses expected information; joint marginal variance
    ## must use the covariance bread (observed for RLB-FC). Older fit objects
    ## have no variance_bread and retain their original expected convention.
    H <- if (!is.null(fit[["variance_bread"]])) fit[["variance_bread"]] else fit$information / n
  } else if (fit$family$family == "poisson") {
    score <- X * (y - p); H <- crossprod(X, X * (w * p)) / n
  } else {
    if (any(p <= 0 | p >= 1)) stop("Log-binomial fit is not interior.")
    score <- X * ((y - p) / (1 - p))
    H <- crossprod(X, X * (w * p * (1 - y) / (1 - p)^2)) / n
  }
  c(target,list(score=score,H=H,covariance=vcov(fit)))
}

maic_delta <- function(terms, v, Z, mu, Vmu, joint, wfit) {
  n <- length(v); v <- v / mean(v)
  point <- maic_target_point(terms,v)
  m0 <- point$risk_C; m1 <- point$risk_A; estimate <- point$estimate
  gradient <- colMeans(v * (terms$dp1 / m1 - terms$dp0 / m0))
  empirical <- v * (terms$p1 / m1 - terms$p0 / m0)
  variance_target <- 0
  diagnostic_message <- ""
  if (!joint) {
    variance_ipd <- drop(crossprod(gradient, terms$covariance %*% gradient))
  } else {
    D <- sweep(Z, 2, mu, "-"); Hg <- sweep(D, 2, colMeans(v * D), "-")
    Ag <- crossprod(D, v * Hg) / n
    # Use Moore-Penrose inverses for both bread blocks, including singular ones.
    ginv_tol <- sqrt(.Machine$double.eps)
    Ai <- MASS::ginv(Ag, tol = ginv_tol); Hi <- MASS::ginv(terms$H, tol = ginv_tol)
    if (any(!is.finite(Ai)) || any(!is.finite(Hi)))
      stop("Nonfinite joint bread inverse.")
    ginv_residual <- function(A, inverse) {
      norm(A %*% inverse %*% A - A, "I") / max(1, norm(A, "I"))
    }
    # Diagnostic only: small reciprocal condition numbers do not discard a SE.
    # Use the existing message field so old/new checkpoint schemas can be mixed.
    diagnostic_message <- sprintf(paste0("Joint bread diagnostics: method=MASS::ginv; tol=%.6g; ",
      "rcond_Ag=%.6g; rcond_H=%.6g; ginv_residual_Ag=%.6g; ginv_residual_H=%.6g."),
      ginv_tol, rcond(Ag), rcond(terms$H), ginv_residual(Ag, Ai), ginv_residual(terms$H, Hi))
    IFg <- -(v * D) %*% t(Ai)
    G <- crossprod(terms$score, wfit * Hg) / n
    IFtheta <- (wfit * terms$score + IFg %*% t(G)) %*% t(Hi)
    gg <- colMeans(Hg * empirical)
    influence <- drop(empirical + IFtheta %*% gradient + IFg %*% gg)
    tg <- drop(t(Hi %*% G %*% Ai) %*% gradient + t(Ai) %*% gg)
    variance_ipd <- sum(influence^2) / n^2
    variance_target <- drop(crossprod(tg, Vmu %*% tg))
  }
  variance <- variance_ipd + variance_target
  se <- if (is.finite(variance) && variance > 0) sqrt(variance) else NA_real_
  list(estimate = estimate, se = se, low = estimate - qnorm(.975) * se,
    up = estimate + qnorm(.975) * se, p = 2 * pnorm(-abs(estimate / se)),
    risk_C = m0, risk_A = m1, variance_ipd = variance_ipd, variance_target = variance_target,
    message = diagnostic_message)
}

maic_fit_status <- function(result, method) {
  if (method == "brmad" && length(result$variables$effect_modifiers))
    return("not_applicable_with_effect_modifiers")
  if (!identical(result$weights$status, "available")) return("weighting_failed")
  fit <- result$fits[[method]]
  if (is.null(fit)) return("fit_failed")
  if (method %in% c("CMH", "CMH_sandwich")) return("available")
  actual <- if (maic_wrapped_glm(method)) fit$fit else fit
  if (is.null(actual)) return("fit_failed")
  if (!isTRUE(actual$converged)) return("fit_not_converged")
  "available"
}

maic_wald_rows <- function(result) {
  v <- result$weights$weight
  vars <- result$variables
  model_maic <- identical(vars$model_weighting, "maic")
  rows <- maic_empty_rows(if (model_maic) "maic" else "unit")
  joint_methods <- c("brm", "brm-FC", "RLP", "RLB",
    "RLB-FC", "RLP-FC", "RLP-BFC")
  prepared <- prepare_brm_score_model(result$data$data.IPD, vars$score_spec,
    vars$effect_modifiers, vars$non_effect_modifiers, vars$model_covariates)
  data <- prepared$data$data; n <- nrow(data)
  for (i in seq_len(nrow(rows))) {
    method <- rows$method[i]
    st <- maic_fit_status(result, method)
    if (st != "available") {
      rows$status[i] <- rows$se_status[i] <- rows$p_status[i] <- st
      next
    }
    if (rows$ci_type[i] == "Exact") next
    t0 <- proc.time()[["elapsed"]]
    fit <- result$fits[[method]]
    ans <- tryCatch({
      if (method %in% c("CMH", "CMH_sandwich")) {
        z <- result$estimates[result$estimates$method == method, , drop = FALSE]
        list(estimate = z$estimate[1], se = z$se[1], low = z$lower[1], up = z$upper[1],
          p = z$p.value[1], status = z$interval_status[1], se_status = z$se_status[1], p_status = z$p_status[1])
      } else {
        wfit <- result$model_weights$weight
        actual <- if (maic_wrapped_glm(method)) fit$fit else fit
        terms <- if (method %in% c("brm", "brm-FC", "brmad"))
          maic_brm_terms(fit, data$y, data$x, wfit) else maic_glm_terms(actual, data, wfit)
        if (maic_wrapped_glm(method)) terms$covariance <- fit$vcov
        joint <- model_maic && method %in% joint_methods
        Z <- mu <- Vmu <- NULL
        if (joint) {
          Z <- as.matrix(data[, vars$weight_covariates, drop = FALSE])
          ad <- result$data$data.AD
          mu <- select_weight_summary(ad$anchored.mean.x, colnames(Z), colnames(Z))
          Vmu <- select_weight_summary(ad$anchored.cov.x, colnames(Z), colnames(Z)) / ad$anchored.n
        }
        maic_delta(terms, v, Z, mu, Vmu, joint, wfit)
      }
    }, error = identity)
    rows$elapsed_seconds[i] <- proc.time()[["elapsed"]] - t0
    if (inherits(ans, "error")) {
      rows$status[i] <- rows$se_status[i] <- rows$p_status[i] <- "standardization_failed"
      rows$message[i] <- conditionMessage(ans)
      # A variance failure must not erase a separately computable GLM target.
      if (method %in% c("LB","LP","RLP","RLB",
                        anchored_glm_fc_method_names())) {
        target <- tryCatch(maic_target_point(maic_glm_target_terms(actual,data),v),error=identity)
        if (!inherits(target,"error")) {
          for (key in names(target)) rows[i,key] <- target[[key]]
          rows$status[i] <- rows$se_status[i] <- rows$p_status[i] <- "se_failed"
        }
      }
    } else {
      for (key in intersect(names(ans), names(rows))) rows[i, key] <- ans[[key]]
      if (!"status" %in% names(ans)) {
        rows$status[i] <- rows$se_status[i] <- rows$p_status[i] <-
          if (is.finite(ans$se) && ans$se > 0) "available" else "se_failed"
      }
    }
    ## Recover the BRM point even if Wald variance failed; match exact's nonlinear target.
    if (method %in% c("brm", "brm-FC", "brmad")) {
      target <- tryCatch(marginal_logrr_target(fit$point.est, fit$va, fit$vb, v, TRUE), error = identity)
      if (!inherits(target, "error")) {
        rows$estimate[i] <- target$value; rows$risk_C[i] <- target$risk.C; rows$risk_A[i] <- target$risk.A
        if (rows$status[i] == "available") {
          rows$low[i] <- target$value - qnorm(.975) * rows$se[i]
          rows$up[i] <- target$value + qnorm(.975) * rows$se[i]
          rows$p[i] <- 2 * pnorm(-abs(target$value / rows$se[i]))
        }
      }
    }
  }
  rows
}

maic_atomic_save <- function(object, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, tmp)
  # POSIX rename replaces the previous checkpoint atomically on the cluster.
  if (.Platform$OS.type == "windows" && file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Cannot finalize result: ", path)
}

maic_estimation_covariates <- function() {
  available <- paste0("z", 1:6)
  parse <- function(value, label) {
    vars <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
    if (!length(vars) || any(!nzchar(vars)) || anyDuplicated(vars) ||
        any(!vars %in% available))
      stop(label, " must list distinct covariates from z1 through z6.")
    vars
  }
  em <- parse(Sys.getenv("MAIC_EST_EM", "z1,z2,z3"), "MAIC_EST_EM")
  non_em_value <- Sys.getenv("MAIC_EST_NON_EM", unset = NA_character_)
  non_em <- if (is.na(non_em_value)) setdiff(available, em) else
    parse(non_em_value, "MAIC_EST_NON_EM")
  if (!length(non_em) || length(intersect(em, non_em)))
    stop("Estimation EM and non-EM sets must be nonempty and disjoint.")
  list(em = em, non_em = non_em)
}

maic_run_seed <- function(code_dir, run_dir, seed, ci_mode = c("full", "smoke"),
                          n = 200L, ess_ratio = .5, fit_control = list()) {
  ci_mode <- match.arg(ci_mode)
  fit_control <- anchored_fit_control(fit_control)
  event <- tolower(Sys.getenv("MAIC_EVENT", "common"))
  if (!event %in% c("common", "rare"))
    stop("MAIC_EVENT must be 'common' or 'rare'.")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 ||
      n != floor(n) || n %% 2 != 0 || n > .Machine$integer.max)
    stop("n must be a positive even integer (total sample size per trial).")
  if (!is.numeric(ess_ratio) || length(ess_ratio) != 1L || !is.finite(ess_ratio) ||
      ess_ratio <= 0 || ess_ratio > 1) stop("ess_ratio must be in (0, 1].")
  n <- as.integer(n)
  weight_mode <- tolower(Sys.getenv("MAIC_WEIGHT_MODE", "all"))
  if (!weight_mode %in% c("all", "em", "nonem"))
    stop("MAIC_WEIGHT_MODE must be 'all', 'em' or 'nonem'.")
  weight_covariates <- paste0("z", switch(weight_mode, all = 1:6, em = 1:3, nonem = 4:6))
  model_weighting <- tolower(Sys.getenv("MAIC_MODEL_WEIGHTING", "unit"))
  if (!model_weighting %in% c("unit", "maic")) stop("MAIC_MODEL_WEIGHTING must be unit or maic.")
  estimation_mode <- tolower(Sys.getenv("MAIC_EST_MODE", "scores"))
  if (!estimation_mode %in% c("scores", "none"))
    stop("MAIC_EST_MODE must be 'scores' or 'none'.")
  if (estimation_mode == "none") {
    estimation_em <- estimation_non_em <- character(0)
    model_covariates <- effect_modifiers <- non_effect_modifiers <- character(0)
    score_spec <- NULL
  } else {
    estimation_covariates <- maic_estimation_covariates()
    estimation_em <- estimation_covariates$em
    estimation_non_em <- estimation_covariates$non_em
    model_covariates <- c("h_E", "h_P")
    effect_modifiers <- "h_E"
    non_effect_modifiers <- "h_P"
    score_spec <- make_brm_score_spec(em_covariates = estimation_em,
      non_em_covariates = estimation_non_em)
  }
  stopifnot(length(seed) == 1L, is.finite(seed), seed >= 1L,
    seed == floor(seed), seed <= .Machine$integer.max - 1000000L)
  controls <- if (ci_mode == "full") list(nsim.coarse = 100L, nsim.confirm = 200L,
    nsim.refine = 200L, nsim.p.value = 500L, time.limit.sec = 10800) else
    list(nsim.coarse = 5L, nsim.confirm = 10L, nsim.refine = 10L, nsim.p.value = 20L,
         search.mult = .5, grid.mult = 25, thres.dicho = .05, bisection.max.step = 8L,
         compute.p.value = FALSE, time.limit.sec = 180)
  dgp_mode <- tolower(Sys.getenv("MAIC_DGP", "legacy"))
  if (!dgp_mode %in% c("legacy", "ess_em")) stop("MAIC_DGP must be legacy or ess_em.")
  if (dgp_mode == "ess_em") {
    design <- webster_rr_scenario_design(event = event, ess_ratio = ess_ratio,
      ess_source = tolower(Sys.getenv("MAIC_ESS_SOURCE", "diffuse")),
      em_level = tolower(Sys.getenv("MAIC_EM_LEVEL", "high")))
  } else {
    design <- webster_rr_design(event = event, ess_ratio = ess_ratio)
  }
  relative_ess <- c(full = rel_ess_binary(design$q, design$e),
    em = rel_ess_binary(design$q[1:3], design$e[1:3]),
    non_em = rel_ess_binary(design$q[4:6], design$e[4:6]))
  ess_source <- if (dgp_mode == "ess_em") design$ess_source else "legacy_diffuse"
  em_level <- if (dgp_mode == "ess_em") design$em_level else "legacy_preset"
  lambda_AC <- unname(design$coefficient_truth["AC", "lambda"])
  lambda_BC <- unname(design$coefficient_truth["BC", "lambda"])
  calibration_ess <- unname(relative_ess[switch(weight_mode, all = "full", em = "em", nonem = "non_em")])
  settings <- list(version = 11L, dgm = design$dgm, dgm_parameters = design$parameters,
    dgp_mode = dgp_mode, dgm_scenario = design$scenario,
    population_q = design$q, population_e = design$e,
    ess_source = ess_source, em_level = em_level, weight_mode = weight_mode,
    theoretical_relative_ess = relative_ess,
    comparison_type = "anchored", method_naming_version = 2L, event = event, n_AC = n, n_BC = n, ess_ratio = ess_ratio,
    model_weighting = model_weighting, cmh_weighting = "maic",
    counterfactual_risk_check = "finite_positive_standardized_means",
    cmh_sandwich_uncertainty = if (model_weighting == "maic") "joint" else "fixed",
    weight_covariates = weight_covariates, op_model = "marginal", score_standardization = "bc_estimated",
    seed = seed, exact_seed = seed + 1000000L, ci_mode = ci_mode, ci_control = controls,
    fit_control = fit_control)
  if (estimation_mode == "none") settings$estimation_mode <- estimation_mode
  if (!identical(estimation_em, paste0("z", 1:3)) ||
      !identical(estimation_non_em, paste0("z", 4:6))) {
    settings$estimation_em_covariates <- estimation_em
    settings$estimation_non_em_covariates <- estimation_non_em
  }
  dir.create(file.path(run_dir, "by_seed"), recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(run_dir, "by_seed", sprintf("seed_%06d", seed))
  if (file.exists(paste0(stem, ".rds"))) {
    old <- readRDS(paste0(stem, ".rds"))
    if (!identical(old$settings, settings)) stop("Existing result has different settings; use a new run directory.")
    if (isTRUE(old$complete)) { cat("Already complete: seed", seed, "\n"); return(invisible(old)) }
  }
  t0 <- proc.time()[["elapsed"]]
  result <- NULL; exact <- list(); warnings <- character(); failure <- NULL
  rows <- maic_empty_rows(model_weighting); truth <- ess <- NA_real_
  source_files <- c("load_anchored.R", unlist(lapply(c("data", "estimation", "summary"),
    function(d) file.path(d, list.files(file.path(code_dir, d), pattern = "\\.[Rr]$", recursive = TRUE)))))
  source_hashes <- tools::md5sum(file.path(code_dir, source_files)); names(source_hashes) <- source_files
  save_checkpoint <- function(complete = FALSE) {
    out <- rows
    out$RR <- exp(out$estimate); out$RR_low <- exp(out$low); out$RR_up <- exp(out$up)
    out$seed <- seed; out$event <- event; out$n_AC <- out$n_BC <- n
    out$dgm <- settings$dgm
    out$ess_ratio <- ess_ratio; out$realized_ess <- ess; out$ci_mode <- ci_mode
    out$realized_ess_ratio <- ess / n
    out$dgp_mode <- dgp_mode; out$ess_source <- ess_source; out$em_level <- em_level
    out$weight_mode <- weight_mode
    out$lambda_AC <- lambda_AC; out$lambda_BC <- lambda_BC; out$lambda_AB <- lambda_AC - lambda_BC
    out$theoretical_ess_ratio_full <- unname(relative_ess["full"])
    out$theoretical_ess_ratio_em <- unname(relative_ess["em"])
    out$theoretical_ess_ratio_nonem <- unname(relative_ess["non_em"])
    out$theoretical_ess_ratio_calibration <- calibration_ess
    out$population_log_RR_AB_in_BC <- unname(design$marginal_truth["BC", "log_RR_AB"])
    out$model_weighting <- model_weighting
    out$weight_covariates <- paste(weight_covariates, collapse = ",")
    out$estimation_mode <- estimation_mode
    out$estimation_em_covariates <- paste(estimation_em, collapse = ",")
    out$estimation_non_em_covariates <- paste(estimation_non_em, collapse = ",")
    out$population_log_RR_AC_in_BC <- ifelse(out$component == "marginal", truth, NA_real_)
    out$seed_complete <- complete
    cf <- NULL
    if (!is.null(result)) {
      cf <- result$estimates[result$estimates$component != "marginal", , drop = FALSE]
      cf$seed <- seed
      cf$event <- event; cf$n_AC <- cf$n_BC <- n; cf$ess_ratio <- ess_ratio
      cf$dgm <- settings$dgm
      cf$dgp_mode <- dgp_mode; cf$ess_source <- ess_source; cf$em_level <- em_level
      cf$weight_mode <- weight_mode
      cf$weight_covariates <- paste(weight_covariates, collapse = ",")
      cf$estimation_mode <- estimation_mode
      cf$estimation_em_covariates <- paste(estimation_em, collapse = ",")
      cf$estimation_non_em_covariates <- paste(estimation_non_em, collapse = ",")
    }
    record <- list(summary = out, coefficients = cf, result = result, exact = exact,
      settings = settings, complete = complete,
      elapsed_seconds = proc.time()[["elapsed"]] - t0,
      warnings = warnings, error = failure, source_md5 = source_hashes,
      job_id = Sys.getenv("SLURM_JOB_ID"), array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID"),
      session_info = sessionInfo())
    # A single atomic RDS checkpoint contains every result for this seed.
    maic_atomic_save(record, paste0(stem, ".rds"))
    invisible(record)
  }
  withCallingHandlers(tryCatch({
    cat("Preparing all methods, n", n, "ESS", ess_ratio, "seed", seed, "mode", ci_mode, "\n")
    truth <- unname(design$marginal_truth["BC", "log_RR_AC"])
    result <- run.anchored(design, n_AC = n, n_BC = n, seed = seed,
      run_exact = FALSE, weight_covariates = settings$weight_covariates,
      model_covariates = model_covariates, effect_modifiers = effect_modifiers,
      non_effect_modifiers = non_effect_modifiers,
      op_model = "marginal", score_standardization = "bc_estimated",
      score_spec = score_spec,
      return_full = TRUE, return_data = TRUE, fit_control = fit_control,
      model_weighting = model_weighting)
    w <- result$weights$weight
    if (all(is.finite(w)) && sum(w^2) > 0) ess <- sum(w)^2 / sum(w^2)
    rows <- maic_wald_rows(result)
    save_checkpoint()
    for (i in which(rows$ci_type == "Exact" & rows$status == "pending")) {
      method <- rows$method[i]; fit <- result$fits[[method]]
      dat <- result$data$data.IPD$data; method_start <- proc.time()[["elapsed"]]
      cat("Starting exact_marginal_target:", method, "seed", seed, "\n"); flush.console()
      ci <- tryCatch(marginal_exact_from_fit(fit, dat$y, dat$x,
        weight = result$model_weights$weight,
        standardization.weight = w, run_exact = TRUE, seed = settings$exact_seed,
        bc_control = controls, label = paste(method, "marginal exact")), error = function(e) {
          z <- make_na_marginal_exact("error"); z$message <- conditionMessage(e); z
        })
      exact[[method]] <- ci
      for (key in c("estimate", "low", "up", "p", "status")) rows[i, key] <- ci[[key]]
      rows$p_status[i] <- ci$p.status; rows$se_status[i] <- "not_provided"
      rows$message[i] <- if (is.null(ci$message)) "" else ci$message
      rows$elapsed_seconds[i] <- proc.time()[["elapsed"]] - method_start
      save_checkpoint()
      cat("Finished", method, "status", ci$status, "seconds", rows$elapsed_seconds[i], "\n"); flush.console()
    }
  }, error = function(e) {
    failure <<- conditionMessage(e)
    pending <- rows$status == "pending"
    rows$status[pending] <<- "pipeline_failed"; rows$message[pending] <<- failure
    rows$p_status[pending] <<- "pipeline_failed"
  }), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning")
  })
  record <- save_checkpoint(complete = is.null(failure) && !any(rows$status == "pending"))
  print(record$summary[c("method", "ci_type", "estimate", "se", "low", "up", "status")])
  if (!is.null(failure)) stop(failure)
  invisible(record)
}

maic_collect <- function(run_dir) {
  spec <- readLines(file.path(run_dir, "submission.txt"), warn = FALSE)
  field <- function(key) sub(paste0("^", key, "="), "", spec[startsWith(spec, paste0(key, "="))])
  expected <- seq.int(as.integer(field("first_seed")), length.out = as.integer(field("nsim")))
  paths <- file.path(run_dir, "by_seed", sprintf("seed_%06d.rds", expected))
  present <- file.exists(paths)
  records <- lapply(paths[present], function(path) {
    record <- readRDS(path)
    if (!is.list(record) || !is.data.frame(record$summary) ||
        !all(c("ci_type", "seed_complete") %in% names(record$summary)))
      stop("Invalid all-method RDS checkpoint: ", path)
    list(summary = record$summary, coefficients = record$coefficients,
         settings = record$settings, complete = isTRUE(record$complete))
  })
  if (length(records) > 1L) {
    first <- records[[1L]]$settings
    first[c("seed", "exact_seed")] <- NULL
    for (record in records[-1L]) {
      current <- record$settings
      current[c("seed", "exact_seed")] <- NULL
      if (!identical(current, first)) stop("Mixed DGP/analysis settings in one result directory.")
    }
  }
  missing <- data.frame(seed = expected[!present])
  write.csv(missing, file.path(run_dir, "missing_seeds.csv"), row.names = FALSE)
  completed <- integer()
  if (length(records)) {
    combined <- do.call(rbind, lapply(records, `[[`, "summary"))
    write.csv(combined, file.path(run_dir, "combined_results.csv"), row.names = FALSE)
    counts <- aggregate(list(n_rows = rep(1L, nrow(combined))), combined[c("method", "ci_type", "status")], sum)
    write.csv(counts, file.path(run_dir, "status_counts.csv"), row.names = FALSE)
    completed <- expected[present][vapply(records, `[[`, logical(1L), "complete")]
    print(counts)
    cf <- Filter(function(x) is.data.frame(x) && nrow(x), lapply(records, `[[`, "coefficients"))
    if (length(cf)) write.csv(do.call(rbind, cf), file.path(run_dir, "combined_coefficients.csv"), row.names = FALSE)
  }
  write.csv(data.frame(seed = setdiff(expected[present], completed)),
    file.path(run_dir, "incomplete_seeds.csv"), row.names = FALSE)
  cat("Recorded:", sum(present), "replicates; complete:", length(completed),
      "; missing:", sum(!present), "; expected:", length(expected), "\n")
  invisible(length(completed))
}
