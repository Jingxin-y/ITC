## MAIC 指数倾斜权重。
## 原代码：compare_RR_anchored.R

## Z 的每列是要平衡的协变量；target_mean 与这些列对应。
get.weight <- function(Z, target_mean, balance_tolerance = 1e-6, maxit = 1000L) {
  Z <- as.matrix(Z)
  if ((!is.numeric(Z) && !is.logical(Z)) || !nrow(Z) || any(!is.finite(Z)))
    stop("Z must be a nonempty finite numeric matrix.")
  if (length(balance_tolerance) != 1L || !is.finite(balance_tolerance) || balance_tolerance <= 0)
    stop("balance_tolerance must be positive and finite.")
  covariates <- colnames(Z)
  q <- ncol(Z)

  if (q == 0) {
    return(list(weight = rep(1, nrow(Z)), ga = setNames(numeric(0), character(0)),
                covariates = character(0), target_mean = numeric(0), status = "available",
                diagnostics = list(converged = TRUE, optimizer_code = 0L,
                  achieved_mean = numeric(0), balance_error = numeric(0),
                  max_scaled_error = 0, tolerance = balance_tolerance)))
  }

  ## 两边都有名称时按矩阵列名对齐；否则按输入顺序对应。
  if (!is.null(covariates) && !is.null(names(target_mean))) {
    target_mean <- target_mean[covariates]
  }
  target_mean <- as.numeric(target_mean)
  if (length(target_mean) != q || any(!is.finite(target_mean))) {
    stop("target_mean must provide one mean for each column of Z.")
  }

  X <- sweep(Z, 2, target_mean, "-")
  scale <- pmax(1, apply(Z, 2, function(z) diff(range(z))))
  fail <- function(reason, message, fit = NULL, achieved = rep(NA_real_, q)) {
    balance <- achieved - target_mean
    diagnostics <- list(converged = FALSE, reason = reason,
      optimizer_code = if (is.null(fit)) NA_integer_ else fit$convergence,
      optimizer_message = if (is.null(fit)) message else fit$message,
      achieved_mean = setNames(achieved, covariates), balance_error = setNames(balance, covariates),
      max_scaled_error = if (all(is.finite(balance))) max(abs(balance) / scale) else NA_real_,
      tolerance = balance_tolerance)
    result <- list(weight = rep(NA_real_, nrow(Z)), ga = setNames(rep(NA_real_, q), covariates),
      covariates = covariates, target_mean = setNames(target_mean, covariates),
      status = "weighting_failed", diagnostics = diagnostics)
    stop(structure(list(message = message, call = NULL, result = result),
                   class = c("maic_calibration_error", "error", "condition")))
  }
  outside <- target_mean < apply(Z, 2, min) - balance_tolerance * scale |
             target_mean > apply(Z, 2, max) + balance_tolerance * scale
  if (any(outside)) fail("target_outside_support", "MAIC target is outside the observed covariate support.",
                         achieved = colMeans(Z))
  Q <- function(alpha, X) log_sum_exp(X %*% alpha)
  Q.grad <- function(alpha, X) {
    eta <- as.vector(X %*% alpha)
    w <- exp(eta - max(eta))
    as.vector(crossprod(X, w / sum(w)))
  }

  fit <- tryCatch(optim(par = rep(0, q), fn = Q, gr = Q.grad, X = X,
               method = "BFGS", control = list(maxit = maxit, reltol = 1e-14)),
               error = function(e) fail("optimizer_error", conditionMessage(e)))

  w <- normalize_exp_weights(X %*% fit$par)
  achieved <- colSums(Z * w) / sum(w)
  balance <- achieved - target_mean
  error <- max(abs(balance) / scale)
  if (fit$convergence != 0L || !is.finite(error) || error > balance_tolerance)
    fail(if (fit$convergence != 0L) "optimizer_not_converged" else "balance_not_achieved",
         sprintf("MAIC calibration failed (optimizer code %s; max scaled balance error %.3g).",
                 fit$convergence, error), fit, achieved)
  list(weight = w, ga = setNames(fit$par, covariates),
       covariates = covariates, target_mean = setNames(target_mean, covariates), status = "available",
       diagnostics = list(converged = TRUE, optimizer_code = fit$convergence,
         achieved_mean = setNames(achieved, covariates), balance_error = setNames(balance, covariates),
         max_scaled_error = error, tolerance = balance_tolerance))
}
