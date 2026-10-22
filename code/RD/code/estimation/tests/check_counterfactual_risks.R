## Untruncated counterfactual predictions must retain computable GLM targets.
## Test real fits, the exporter, point recovery and essential numerical checks.
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
source("estimation/cluster_all_methods.R", encoding = "UTF-8")

close_to <- function(a, b, tol = 1e-6)
  stopifnot(isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE)))
probability_methods <- c("LB", "RLB", "LB-FC", "RLB-FC", "LP-BFC", "RLP-BFC")
poisson_methods <- c("LP", "RLP", "LP-FC", "RLP-FC")
glm_methods <- c(probability_methods, poisson_methods)
actual_fit <- function(result, method) {
  fit <- result$fits[[method]]
  if (maic_wrapped_glm(method)) fit$fit else fit
}

make_result <- function(valid, mode, fit_control = list()) {
  cells <- if (valid) data.frame(x = c(0, 0, 1, 1), z1 = c(0, 1, 0, 1),
    events = c(10, 20, 30, 60)) else
    data.frame(x = c(0, 0, 1), z1 = c(0, 1, 0), events = c(20, 80, 60))
  dat <- do.call(rbind, lapply(seq_len(nrow(cells)), function(k)
    data.frame(id = 1:100, x = cells$x[k], z1 = cells$z1[k],
      y = as.integer(1:100 <= cells$events[k]))))
  ipd <- list(data = dat)
  ## A 1% target prevalence keeps BOTH averaged LB risks below 1 even when
  ## 100 individual A risks exceed 1. Point computability is a separate issue.
  ad <- get.aggre(list(data = data.frame(z1 = rep(c(0, 1), c(99, 1)))))
  z <- suppressWarnings(get.compare.anchored(ipd, ad, run_exact = FALSE,
    return_full = TRUE, score_standardization = "fixed", weight_covariates = "z1",
    model_covariates = "z1", effect_modifiers = character(),
    non_effect_modifiers = "z1", model_weighting = mode, fit_control = fit_control))
  z$data <- list(data.IPD = ipd, data.AD = ad)
  z
}

check_recovery <- function(valid, invalid) {
  real_delta <- maic_delta
  on.exit(assign("maic_delta", real_delta, envir = .GlobalEnv))
  assign("maic_delta", function(...) stop("Injected variance failure."), envir = .GlobalEnv)
  valid_rows <- maic_wald_rows(valid)
  good <- valid_rows$method %in% glm_methods
  stopifnot(all(valid_rows$status[good] == "se_failed"),
    all(is.finite(valid_rows$estimate[good])), all(is.na(valid_rows$se[good])))
  invalid_rows <- maic_wald_rows(invalid)
  computable <- invalid_rows$method %in% glm_methods
  stopifnot(all(invalid_rows$status[computable] == "se_failed"),
    all(is.finite(invalid_rows$estimate[computable])),
    all(is.na(invalid_rows$se[computable])))
}

for (mode in c("unit", "maic")) {
  invalid <- make_result(FALSE, mode)
  dat <- invalid$data$data.IPD$data
  rows <- maic_wald_rows(invalid)
  for (method in glm_methods) {
    fit <- actual_fit(invalid, method)
    stopifnot(isTRUE(fit$converged), all(fitted(fit) > 0 & fitted(fit) < 1))
    ## Conditional coefficients still describe the converged factual fit.
    coefficients <- invalid$estimates[invalid$estimates$method == method &
      invalid$estimates$component != "marginal", , drop = FALSE]
    stopifnot(nrow(coefficients) > 0, all(coefficients$status == "available"))
    target <- maic_glm_target_terms(fit, dat)
    stopifnot(sum(target$p1 > 1) == 100L)
    d0 <- d1 <- dat; d0$x <- 0; d1$x <- 1
    p0 <- predict(fit, d0, type = "response"); p1 <- predict(fit, d1, type = "response")
    close_to(target$p0, p0); close_to(target$p1, p1)
    means <- c(weighted.mean(p0, invalid$weights$weight),
      weighted.mean(p1, invalid$weights$weight))
    row <- rows[rows$method == method & rows$ci_type == "Wald", ]
    close_to(c(row$risk_C, row$risk_A), means)
    close_to(row$estimate, log(means[2] / means[1]))
  }
  computable <- rows$method %in% glm_methods
  stopifnot(all(rows$status[computable] == "available"),
    all(is.finite(rows$estimate[computable])), all(is.finite(rows$se[computable])),
    all(is.finite(rows$low[computable])), all(is.finite(rows$up[computable])))

  fit <- actual_fit(invalid, "LB")
  d0 <- d1 <- dat; d0$x <- 0; d1$x <- 1
  p0 <- predict(fit, d0, type = "response"); p1 <- predict(fit, d1, type = "response")
  close_to(max(p1), 2.4)
  means <- c(weighted.mean(p0, invalid$weights$weight), weighted.mean(p1, invalid$weights$weight))
  close_to(means, c(.206, .618))
  stopifnot(all(means > 0 & means < 1))
  ## A standardized working mean above 1 also remains computable.
  point <- maic_target_point(maic_glm_target_terms(fit, dat), rep(1, nrow(dat)))
  close_to(c(point$risk_C, point$risk_A), c(.4, 1.2))
  close_to(point$estimate, log(3))

  valid <- make_result(TRUE, mode)
  valid_dat <- valid$data$data.IPD$data
  valid_rows <- maic_wald_rows(valid)
  stopifnot(all(valid_rows$status[valid_rows$method %in% glm_methods] == "available"))
  for (method in glm_methods) {
    fit <- actual_fit(valid, method)
    target <- maic_glm_target_terms(fit, valid_dat)
    d0 <- d1 <- valid_dat; d0$x <- 0; d1$x <- 1
    p0 <- predict(fit, d0, type = "response"); p1 <- predict(fit, d1, type = "response")
    close_to(target$p0, p0); close_to(target$p1, p1)
    point <- c(risk_C = weighted.mean(p0, valid$weights$weight),
      risk_A = weighted.mean(p1, valid$weights$weight))
    row <- valid_rows[valid_rows$method == method & valid_rows$ci_type == "Wald", ]
    close_to(c(row$risk_C, row$risk_A), point)
    close_to(row$estimate, log(point["risk_A"] / point["risk_C"]))
  }
  check_recovery(valid, invalid)
  ## Probability 1 is no longer a reason to discard a counterfactual target.
  boundary_fit <- actual_fit(valid, "LB")
  boundary_fit$coefficients[] <- 0
  point <- maic_target_point(maic_glm_target_terms(boundary_fit, valid_dat), valid$weights$weight)
  close_to(c(point$risk_C, point$risk_A, point$estimate), c(1, 1, 0))
  ## Nonfinite coefficients and noncomputable standardized means still fail.
  boundary_fit$coefficients[1] <- Inf
  err <- tryCatch(maic_glm_target_terms(boundary_fit, valid_dat), error = identity)
  stopifnot(inherits(err, "error"), grepl("Nonfinite", conditionMessage(err), fixed = TRUE))
  for (intercept in c(-1000, 1000)) {
    boundary_fit <- actual_fit(valid, "LB")
    boundary_fit$coefficients[] <- 0
    boundary_fit$coefficients["(Intercept)"] <- intercept
    err <- tryCatch(maic_target_point(maic_glm_target_terms(boundary_fit, valid_dat),
      valid$weights$weight), error = identity)
    stopifnot(inherits(err, "error"),
      grepl("Invalid standardized risk", conditionMessage(err), fixed = TRUE))
  }
  ## Removing the prediction filter must not relax the Bernoulli FC fitter.
  fit <- actual_fit(valid, "LB-FC")
  theta <- coef(fit); theta[] <- 0
  err <- tryCatch(glm_fc_score(theta, model.matrix(fit), valid_dat$y,
    valid$model_weights$weight, "binomial", "bernoulli"), error = identity)
  stopifnot(inherits(err, "error"),
    grepl("FC iterate is outside", conditionMessage(err), fixed = TRUE))
  failed <- make_result(FALSE, mode, list(fc_max_step = 1L, fc_thres = 1e-14))
  nonconverged <- anchored_glm_fc_method_names()[vapply(anchored_glm_fc_method_names(),
    function(method) identical(maic_fit_status(failed, method), "fit_not_converged"), logical(1))]
  stopifnot(length(nonconverged) > 0L)
  failed_rows <- maic_wald_rows(failed)
  rejected <- failed_rows$method %in% nonconverged
  stopifnot(all(failed_rows$status[rejected] == "fit_not_converged"),
    all(is.na(failed_rows$estimate[rejected])))
  cat("PASS:", mode, "unfiltered targets, numerical guards, factual FC domain and point recovery\n")
}
