## Preserve expected parameters even when a fit fails; track each output separately.
comparison_schema <- function(data, model_covariates, effect_modifiers, op_design) {
  pa <- ncol(op_design$va); pb <- ncol(op_design$vb)
  brm <- data.frame(component = c(rep("alpha", pa), rep("beta", pb)),
    parameter = paste(c(rep("alpha", pa), rep("beta", pb)), c(seq_len(pa), seq_len(pb))),
    term = c(colnames(op_design$va), colnames(op_design$vb)))
  terms <- colnames(model.matrix(weighted_glm_formula(data, model_covariates, effect_modifiers), data))
  glm <- data.frame(component = ifelse(terms == "x" | grepl("^x:", terms), "treatment", "baseline"),
                    parameter = terms, term = terms)
  cmh <- data.frame(component = "common_effect", parameter = "log_RR", term = "common log RR")
  list(brm = brm, glm = glm, robust_glm = glm, cmh = cmh)
}

bc_parameter_status <- function(bc, pa, run_exact = TRUE) {
  out <- data.frame(interval_status = rep(if (run_exact) "BC_failed" else "BC_not_run", pa),
                    p_status = rep(if (run_exact) "BC_failed" else "BC_not_run", pa))
  if (!run_exact || is.null(bc)) return(out)
  failed_fits <- function(d) any(vapply(c("failed.nuisance.fits", "optimizer.failures"),
    function(k) !is.null(d[[k]]) && (!is.finite(d[[k]]) || d[[k]] > 0), logical(1)))
  endpoint_status <- function(d) {
    if (is.null(d) || is.null(d$convergence)) return("BC_unverified")
    if ((!is.null(d$candidate.converged) && any(!d$candidate.converged)) ||
        failed_fits(d$diagnostics)) return("BC_fit_failed")
    if (!isTRUE(d$convergence)) {
      if (identical(d$reason, "confirmed rejected point not found in search range"))
        return("BC_boundary_truncated")
      return("BC_search_failed")
    }
    "available"
  }
  at <- function(x, j) if (length(x) >= j) x[[j]] else NULL
  for (j in seq_len(pa)) {
    statuses <- c(endpoint_status(at(bc$diagnostics$lower, j)),
                  endpoint_status(at(bc$diagnostics$upper, j)))
    accepted <- statuses %in% c("available","BC_boundary_truncated")
    out$interval_status[j] <- if (!all(accepted)) statuses[which(!accepted)[1]] else
      if (any(statuses == "BC_boundary_truncated")) "BC_boundary_truncated" else "available"
    bounds <- c(at(bc$low, j), at(bc$up, j))
    if (out$interval_status[j] %in% c("available","BC_boundary_truncated") &&
        (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1] > bounds[2]))
      out$interval_status[j] <- "BC_invalid_interval"
    pdiag <- at(bc$diagnostics$p.value, j)
    null <- at(bc$diagnostics$observed.null, j)
    pv <- at(bc$p, j)
    out$p_status[j] <- if (is.null(pdiag) || is.null(null$converged)) "BC_unverified" else
      if (!isTRUE(null$converged) || failed_fits(pdiag$diagnostics)) "BC_fit_failed" else
      if (length(pv) != 1L || !is.finite(pv) || pv < 0 || pv > 1) "BC_invalid_p" else "available"
  }
  out
}

finish_anchored_comparison <- function(fits, intervals, weights, variables, schema,
                                       run_exact, return_full, model_weights = weights) {
  kind_by_method <- c(brm = "brm", LB = "glm", LP = "glm", RLP = "robust_glm",
    CMH = "cmh", CMH_sandwich = "cmh", brmad = "brm", "brm-BC" = "brm",
    "brm-FC-BC" = "brm", "brm_b-BC" = "brm", "brm-FC" = "brm", RLB = "robust_glm",
    setNames(rep("robust_glm", length(anchored_glm_fc_method_names())),
             anchored_glm_fc_method_names()))
  kinds <- unname(kind_by_method[names(fits)])
  if (anyNA(kinds) || anyDuplicated(names(fits))) stop("Invalid anchored method names.")
  interval_by_method <- c("brm-BC" = "brm.bc", "brm-FC-BC" = "firth.bc", "brm_b-BC" = "bayes.bc")
  use_bc <- names(fits) %in% names(interval_by_method)
  bc_for <- function(i) if (use_bc[i]) intervals[[interval_by_method[[names(fits)[i]]]]] else NULL
  marginal_bc <- identical(variables$bc_target,"marginal")
  rows <- lapply(seq_along(fits), function(i) {
    method <- names(fits)[i]
    status <- if (length(variables$effect_modifiers) && method %in% c("brmad", "brm_b-BC"))
      "not_applicable_with_effect_modifiers" else
      if (model_weights$status != "available") "weighting_failed" else "available"
    ## Unit-weight coefficients remain valid if only target calibration fails.
    if ((kinds[i] == "cmh" || (marginal_bc && use_bc[i])) &&
        status == "available" && weights$status != "available")
      status <- "weighting_failed"
    if (marginal_bc && use_bc[i]) return(comparison_marginal_rows(
      method,fits[[i]],bc_for(i),weights$weight,run_exact,status))
    comparison_parameter_rows(method, fits[[i]], kinds[i], bc = bc_for(i),
      use_bc = use_bc[i], run_exact = run_exact, status = status, schema = schema[[kinds[i]]])
  })
  estimates <- do.call(rbind, rows)
  nonconverged <- vapply(fits,function(f) isFALSE(f[["converged"]]) ||
    (is.list(f[["fit"]]) && isFALSE(f[["fit"]][["converged"]])),logical(1))
  candidates <- lapply(which(nonconverged),function(i) {
    if (marginal_bc && use_bc[i]) return(comparison_marginal_rows(
      names(fits)[i],fits[[i]],bc_for(i),weights$weight,run_exact,
      retain_nonconverged=TRUE))
    comparison_parameter_rows(names(fits)[i],fits[[i]],kinds[i],schema=schema[[kinds[i]]],
                              retain_nonconverged=TRUE)
  })
  nonconverged_estimates <- if (length(candidates)) do.call(rbind,candidates) else estimates[FALSE,,drop=FALSE]
  if (return_full) return(list(estimates = estimates, fits = fits, bc = intervals,
    marginal_estimates = estimates[estimates$component=="marginal",,drop=FALSE],
    weights = weights, model_weights = model_weights, variables = variables,
    nonconverged_estimates = nonconverged_estimates, nonconverged_fits = fits[nonconverged],
    parameter_schema = unique(estimates[c("method", "component", "parameter", "term")])))
  result <- vapply(seq_along(rows), function(i) {
    d <- rows[[i]]
    selected <- if (marginal_bc && use_bc[i]) 1L else if (kinds[i] == "brm") which(d$component == "alpha")[1] else
      if (kinds[i] %in% c("glm", "robust_glm")) which(d$term == "x")[1] else 1L
    as.numeric(d[selected, c("estimate", "se", "lower", "upper", "p.value")])
  }, numeric(5))
  dimnames(result) <- list(c("point.est", "se.est", "con.low", "con.up", "p.value"), names(fits))
  attr(result, "bc") <- intervals
  attr(result, "estimates") <- estimates
  attr(result, "marginal_estimates") <- estimates[estimates$component=="marginal",,drop=FALSE]
  attr(result, "weights") <- weights
  attr(result, "model_weights") <- model_weights
  attr(result, "nonconverged_estimates") <- nonconverged_estimates
  attr(result, "nonconverged_fits") <- fits[nonconverged]
  result
}

comparison_parameter_rows <- function(method, fit, kind, bc = NULL,
                                      use_bc = FALSE, run_exact = TRUE,
                                      status = "available", schema = NULL,
                                      retain_nonconverged = FALSE) {
  if (is.null(schema)) {
    if (kind == "brm" && !is.null(fit$coefficients)) {
      schema <- data.frame(component = c(rep("alpha", ncol(fit$va)), rep("beta", ncol(fit$vb))),
        parameter = rownames(fit$coefficients), term = c(colnames(fit$va), colnames(fit$vb)))
    } else if (kind == "brm" && isTRUE(fit$bayes_adjusted)) {
      schema <- data.frame(component="alpha",parameter="alpha 1",term="(Intercept)")
    } else if (kind %in% c("glm", "robust_glm") && !is.null(fit)) {
      terms <- names(if (kind == "glm") coef(fit) else fit$point.est)
      schema <- data.frame(component = ifelse(terms == "x" | grepl("^x:", terms),
        "treatment", "baseline"), parameter = terms, term = terms)
    } else if (kind == "cmh") {
      schema <- data.frame(component = "common_effect", parameter = "log_RR", term = "common log RR")
    } else schema <- data.frame(component = NA_character_, parameter = NA_character_, term = NA_character_)
  }
  nr <- nrow(schema)
  out <- cbind(data.frame(method = rep(method, nr)), schema,
    data.frame(estimate = rep(NA_real_, nr), se = NA_real_, lower = NA_real_, upper = NA_real_,
      p.value = NA_real_, interval = if (use_bc) "BC" else "Wald", status = status,
      se_status = status, interval_status = status, p_status = status))
  if (status != "available") return(out)
  missing_fit <- is.null(fit) || (kind == "brm" && is.null(fit$coefficients) && !isTRUE(fit$bayes_adjusted)) ||
    (kind == "robust_glm" && is.null(fit$point.est))
  if (missing_fit) {
    out[c("status", "se_status", "interval_status", "p_status")] <- "fit_failed"
    return(out)
  }
  converged <- if (kind == "glm") isTRUE(fit$converged) else
    if (kind == "brm") isTRUE(fit$converged) else
    if (kind == "robust_glm") isTRUE(fit$fit$converged) else TRUE
  out$status <- if (converged) "available" else "fit_not_converged"
  if (kind == "brm") {
    cf <- fit$coefficients
    if (!is.null(cf)) {
      idx <- match(out$parameter, rownames(cf))
      out[c("estimate", "se", "lower", "upper", "p.value")] <- unname(cf[idx, 1:5, drop = FALSE])
    }
    if (isTRUE(fit$bayes_adjusted)) {
      alpha <- which(out$component == "alpha")
      if (is.null(cf) && length(alpha) == 1L)
        out[alpha,c("estimate", "se", "lower", "upper", "p.value")] <-
          list(fit$point.est[1],fit$se.est[1],fit$conf.lower[1],fit$conf.upper[1],fit$p.value[1])
      out$status[alpha] <- "available"
      out$interval[alpha] <- "Beta posterior"
    }
  } else if (kind %in% c("glm", "robust_glm")) {
    point <- if (kind == "glm") coef(fit) else fit$point.est
    V <- if (kind == "glm") vcov(fit) else fit$vcov
    se <- sqrt(diag(V)); idx <- match(out$term, names(point))
    critical <- if (kind == "glm") qnorm(.975) else 1.96
    out$estimate <- unname(point[idx]); out$se <- unname(se[idx])
    out$lower <- out$estimate - critical * out$se
    out$upper <- out$estimate + critical * out$se
    out$p.value <- 2 * pnorm(-abs(out$estimate / out$se))
  } else out[c("estimate", "se", "lower", "upper", "p.value")] <-
    list(fit$est, fit$se, fit$low, fit$up, fit$p)
  out$status[out$status == "available" & !is.finite(out$estimate)] <- "estimate_nonfinite"
  out$se_status <- ifelse(is.finite(out$se) & out$se > 0, "available", "se_failed")
  out$interval_status <- ifelse(out$se_status == "available" & is.finite(out$lower) &
    is.finite(out$upper) & out$lower <= out$upper, "available", "interval_failed")
  out$p_status <- ifelse(out$se_status == "available" & is.finite(out$p.value) &
    out$p.value >= 0 & out$p.value <= 1, "available", "p_failed")
  if (use_bc) {
    alpha <- which(out$component == "alpha")
    bs <- bc_parameter_status(bc, length(alpha), run_exact)
    out$interval[alpha] <- if (run_exact) "BC" else "BC_not_run"
    out$interval_status[alpha] <- bs$interval_status; out$p_status[alpha] <- bs$p_status
    out$lower[alpha] <- out$upper[alpha] <- out$p.value[alpha] <- NA_real_
    for (j in seq_along(alpha)) {
      if (bs$interval_status[j] %in% c("available","BC_boundary_truncated")) {
        out$lower[alpha[j]] <- bc$low[j]; out$upper[alpha[j]] <- bc$up[j]
      }
      if (bs$p_status[j] == "available") out$p.value[alpha[j]] <- bc$p[j]
    }
  }
  failed <- out$status != "available"
  for (field in c("se_status", "interval_status", "p_status")) out[[field]][failed] <- out$status[failed]
  if (!retain_nonconverged) {
    out$estimate[failed] <- NA_real_
    out$se[out$se_status != "available"] <- NA_real_
    out$lower[!out$interval_status %in% c("available","BC_boundary_truncated")] <- NA_real_
    out$upper[!out$interval_status %in% c("available","BC_boundary_truncated")] <- NA_real_
    out$p.value[out$p_status != "available"] <- NA_real_
  }
  rownames(out) <- NULL
  out
}
