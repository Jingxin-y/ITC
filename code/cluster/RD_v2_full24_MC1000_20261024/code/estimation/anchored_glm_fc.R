## Mean-bias reduction for log-binomial and log-Poisson estimating equations.
## Analytic weights describe independent subjects, not replicated observations.
## See GLM_FC_METHODS.md for the model assumptions and variance conventions.
anchored_glm_fc_method_names <- function() c(
  "LB-FC", "RLB-FC", "LP-FC", "RLP-FC", "LP-BFC", "RLP-BFC")

glm_fc_score <- function(theta, X, y, w, family = c("binomial", "poisson"),
                         bias_model = c("bernoulli", "poisson")) {
  family <- match.arg(family); bias_model <- match.arg(bias_model)
  if (family == "binomial" && bias_model != "bernoulli")
    stop("Log-binomial FC requires the Bernoulli bias model.")
  mu <- exp(drop(X %*% theta))
  if (any(!is.finite(mu)) || any(mu <= 1e-12) ||
      (bias_model == "bernoulli" && any(mu >= 1 - 1e-10)))
    stop("FC iterate is outside the untruncated mean domain.")
  if (family == "binomial") {
    curvature <- mu / (1 - mu)
    score_rows <- X * ((y - mu) / (1 - mu))
    meat_scale <- curvature
  } else {
    curvature <- mu
    score_rows <- X * (y - mu)
    meat_scale <- if (bias_model == "poisson") mu else mu * (1 - mu)
  }
  information <- crossprod(X, X * (w * curvature))
  information.inv <- solve(information)
  score.variance <- crossprod(X, X * (w^2 * meat_scale))
  covariance <- information.inv %*% score.variance %*% t(information.inv)
  leverage.v <- rowSums((X %*% covariance) * X)
  adjustment_rows <- .5 * w * curvature * leverage.v
  if (family == "binomial") {
    leverage.i <- rowSums((X %*% information.inv) * X)
    adjustment_rows <- adjustment_rows +
      w * mu^2 / (1 - mu)^2 * (leverage.v - w * leverage.i)
  }
  score <- colSums(w * score_rows)
  adjustment <- drop(crossprod(X, adjustment_rows))
  list(mu = mu, score = score, adjustment = adjustment,
    adjusted.score = score + adjustment, score_rows = score_rows,
    information = information, information.inv = information.inv,
    score.variance = score.variance, covariance = (covariance + t(covariance)) / 2)
}

solve_glm_fc <- function(start, X, y, w, family, bias_model, maxit, tolerance) {
  theta <- start; evaluations <- 0L
  evaluate <- function(b) {
    evaluations <<- evaluations + 1L
    z <- glm_fc_score(b, X, y, w, family, bias_model)
    if (any(!is.finite(c(z$adjusted.score, z$covariance)))) stop("Nonfinite FC equation.")
    z
  }
  candidate <- function(b) tryCatch(evaluate(b), error = function(e) NULL)
  current <- evaluate(theta)
  scale <- pmax(1, sqrt(diag(current$information)))
  residual <- function(z) z$adjusted.score / scale
  r <- residual(current); step <- 0L; reason <- "maximum iterations reached"
  while (max(abs(r)) > tolerance && step < maxit) {
    step <- step + 1L
    J <- matrix(NA_real_, length(theta), length(theta))
    for (j in seq_along(theta)) {
      h <- 1e-5 * max(1, abs(theta[j])); plus <- minus <- theta
      plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
      fp <- candidate(plus); fm <- candidate(minus)
      if (!is.null(fp) && !is.null(fm)) J[, j] <- (residual(fp) - residual(fm)) / (2*h) else
        if (!is.null(fp)) J[, j] <- (residual(fp) - r) / h else
          if (!is.null(fm)) J[, j] <- (r - residual(fm)) / h
    }
    if (any(!is.finite(J))) { reason <- "could not differentiate adjusted score"; break }
    direction <- tryCatch(drop(solve(J, -r)), error = function(e) -drop(crossprod(J, r)))
    if (any(!is.finite(direction)) || max(abs(direction)) == 0) {
      reason <- "no finite descent direction"; break
    }
    direction <- direction / max(1, max(abs(direction)) / 2)
    loss <- sum(r^2); accepted <- FALSE
    for (halving in 0:24) {
      fraction <- 2^(-halving); next.theta <- theta + fraction * direction
      z <- candidate(next.theta)
      if (is.null(z)) next
      nr <- residual(z)
      if (sum(nr^2) < loss * (1 - 1e-4*fraction) || max(abs(nr)) <= tolerance) {
        theta <- next.theta; current <- z; r <- nr; accepted <- TRUE; break
      }
    }
    if (!accepted) { reason <- "line search did not reduce adjusted-score residual"; break }
  }
  converged <- max(abs(r)) <= tolerance
  if (converged) reason <- "adjusted-score tolerance reached"
  list(par = theta, converged = converged, terms = current,
    diagnostics = list(solver = "joint_damped_newton", reason = reason,
      iterations = step, evaluations = evaluations, max_scaled_score = max(abs(r)),
      score_tolerance = tolerance, adjusted_score = current$adjusted.score,
      score_scale = scale, bias_model = bias_model,
      weights_normalized_to_mean_one = TRUE))
}

## Standard GLM accessors for the reduced-bias fit; no refitting at plain MLE.
coef.glm_fc_fit <- function(object, ...) object$coefficients
vcov.glm_fc_fit <- function(object, ...) object$covariance
fitted.glm_fc_fit <- function(object, ...) object$fitted.values
model.matrix.glm_fc_fit <- function(object, ...) object$x
terms.glm_fc_fit <- function(x, ...) x$terms
predict.glm_fc_fit <- function(object, newdata = NULL, type = c("link", "response"), ...) {
  type <- match.arg(type)
  X <- if (is.null(newdata)) object$x else model.matrix(delete.response(object$terms),
    newdata, contrasts.arg = object$contrasts, xlev = object$xlevels)
  eta <- drop(X %*% object$coefficients)
  if (type == "response") exp(eta) else eta
}

fit_logglm_fc <- function(data, w = NULL, family = c("binomial", "poisson"),
                          bias_model = NULL, covariance = c("model", "sandwich"),
                          model_covariates = NULL, effect_modifiers = NULL,
                          gamma = NULL, target_mean = NULL, target_var = NULL,
                          m_AD = NULL, target_cov = NULL, target_mean_vcov = NULL,
                          account_for_gamma = FALSE, account_for_target_mean = FALSE,
                          weight_covariates = NULL, maxit = 300L, tolerance = 1e-6,
                          start = NULL) {
  family <- match.arg(family); covariance <- match.arg(covariance)
  if (is.null(bias_model)) bias_model <- if (family == "binomial") "bernoulli" else "poisson"
  bias_model <- match.arg(bias_model, c("bernoulli", "poisson"))
  n <- nrow(data)
  if (is.null(w)) w <- rep(1, n)
  if (!n || length(w) != n || any(!is.finite(w)) || any(w < 0) || !any(w > 0))
    stop("FC requires finite nonnegative weights and positive total weight.")
  w <- w / max(w); w <- w / mean(w)
  formula <- weighted_glm_formula(data, model_covariates, effect_modifiers)
  frame <- model.frame(formula, data, na.action = na.fail)
  X <- model.matrix(formula, frame); y <- model.response(frame)
  if (any(!is.finite(c(X, y))) || any(!y %in% c(0, 1)) ||
      any(!is.finite(data$x)) || any(!data$x %in% c(0,1)))
    stop("These anchored FC methods require complete finite binary observations.")
  if (qr(X[w > 0, , drop = FALSE])$rank != ncol(X))
    stop("FC design is not identifiable on positive-weight observations.")
  if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1 || maxit != floor(maxit) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("FC requires a positive integer maxit and positive tolerance.")
  if (account_for_target_mean && !account_for_gamma)
    stop("Target-mean uncertainty requires account_for_gamma=TRUE.")
  if (account_for_gamma && covariance != "sandwich")
    stop("Calibration propagation is available in the sandwich FC versions.")
  if (is.null(start)) {
    start <- setNames(rep(0, ncol(X)), colnames(X))
    for (arm in 0:1) {
      ix <- data$x == arm & w > 0
      if (!any(ix)) stop("Both treatment arms require positive fitting weight.")
      p <- min(.8, max(1e-6, (sum(w[ix]*y[ix])+.5)/(sum(w[ix])+1)))
      if (arm == 0) start["(Intercept)"] <- log(p) else start["x"] <- log(p)-start["(Intercept)"]
    }
    ## A finite likelihood iterate is only an initial value. Root convergence
    ## is still assessed against the adjusted equation, including boundary cases.
    initial <- tryCatch(suppressWarnings(glm(formula, data = frame, weights = w,
      family = if (family == "binomial") binomial("log") else poisson("log"),
      start = start, control = glm.control(maxit = 100L))), error = function(e) NULL)
    if (!is.null(initial) && all(is.finite(coef(initial))) &&
        max(abs(coef(initial))) < 10 &&
        min(exp(drop(X %*% coef(initial)))) > 1e-5) {
      start <- coef(initial)
      max.mu <- max(exp(drop(X %*% start)))
      if (bias_model == "bernoulli" && max.mu >= .95)
        start["(Intercept)"] <- start["(Intercept)"] + log(.95/max.mu)
    }
  } else {
    if (!is.null(names(start))) start <- start[colnames(X)]
    if (length(start) != ncol(X) || any(!is.finite(start))) stop("Invalid FC start.")
  }
  sol <- solve_glm_fc(start, X, y, w, family, bias_model, as.integer(maxit), tolerance)
  theta <- setNames(sol$par, colnames(X)); tt <- sol$terms
  sol$diagnostics$min_fitted_mean <- min(tt$mu)
  sol$diagnostics$max_fitted_mean <- max(tt$mu)
  sol$diagnostics$n_fitted_means_above_one <- sum(tt$mu > 1)
  V <- Vipd <- matrix(NA_real_, ncol(X), ncol(X)); Vtarget <- matrix(0, ncol(X), ncol(X))
  ## Bread of the original score, averaged over subjects. Keep the expected
  ## information in tt unchanged: it still defines the FC point adjustment.
  observed_bread <- family == "binomial" && covariance == "sandwich"
  variance_bread <- if (observed_bread)
    crossprod(X, X * (w * tt$mu * (1 - y) / (1 - tt$mu)^2)) / n else
    tt$information / n
  variance_bread_type <- if (observed_bread) "observed" else "expected"
  influence <- NULL; variance_status <- if (sol$converged) "available" else "fit_not_converged"
  if (sol$converged) {
    variance <- tryCatch({
      if (covariance == "model") Vipd <- tt$covariance else {
        ## First-order covariance: the O(1) score adjustment does not change
        ## the leading influence function of the unadjusted estimating equation.
        ## RLB-FC follows RLB's observed derivative and Moore-Penrose inverse.
        ## Log-Poisson keeps its existing (observed = expected) bread.
        Hi <- if (observed_bread) MASS::ginv(variance_bread) else tt$information.inv * n
        influence <- (w * tt$score_rows) %*% t(Hi)
        if (account_for_gamma) {
          cols <- weight_covariate_names(data, weight_covariates, gamma)
          Z <- covariate_matrix(data, cols); q <- ncol(Z)
          if (q) {
            mu <- select_weight_summary(target_mean, cols, covariate_names(data))
            if (is.null(mu) || length(mu) != q || any(!is.finite(mu)))
              stop("FC joint covariance needs target means for calibration variables.")
            if (!is.null(gamma)) {
              ga <- if (is.null(names(gamma))) gamma else gamma[cols]
              if (length(ga) != q || any(!is.finite(ga))) stop("Invalid FC calibration gamma.")
              reconstructed <- normalize_exp_weights(sweep(Z,2,mu,"-") %*% ga)
              if (max(abs(reconstructed-w)) > 1e-6*max(1,w))
                stop("FC joint covariance requires gamma matching the supplied weights.")
            }
            D <- sweep(Z, 2, mu, "-"); Hg <- sweep(D, 2, colSums(w*D)/sum(w), "-")
            Agi <- MASS::ginv(crossprod(D, w*Hg)/n)
            IFg <- -(w*D) %*% t(Agi)
            G <- crossprod(tt$score_rows, w*Hg)/n
            influence <- (w*tt$score_rows + IFg %*% t(G)) %*% t(Hi)
            if (account_for_target_mean) {
              Vm <- make_target_mean_vcov(
                select_weight_summary(target_var, cols, covariate_names(data)), m_AD,
                select_weight_summary(target_cov, cols, covariate_names(data)),
                select_weight_summary(target_mean_vcov, cols, covariate_names(data)))
              if (!all(dim(Vm) == c(q, q))) stop("Invalid target-mean covariance dimension.")
              sensitivity <- Hi %*% G %*% Agi * mean(w)
              Vtarget <- sensitivity %*% Vm %*% t(sensitivity)
            }
          }
        }
        Vipd <- crossprod(influence)/n^2
      }
      V <- (Vipd + Vtarget + t(Vipd + Vtarget))/2
      if (any(!is.finite(V)) || any(diag(V) <= 0)) stop("Invalid FC covariance.")
      TRUE
    }, error = function(e) e)
    if (inherits(variance, "error")) {
      variance_status <- "se_failed"; V[,] <- NA_real_
      warning("FC covariance failed: ", conditionMessage(variance), call. = FALSE)
    }
  } else warning("GLM FC did not converge; last iterate retained separately.", call. = FALSE)
  dimnames(V) <- dimnames(Vipd) <- dimnames(Vtarget) <- list(colnames(X), colnames(X))
  fit <- structure(list(coefficients = theta, fitted.values = tt$mu,
    linear.predictors = drop(X %*% theta), prior.weights = w, y = y, x = X,
    terms = terms(formula, data = data), formula = formula, model = frame,
    contrasts = attr(X, "contrasts"), xlevels = stats::.getXlevels(terms(formula), frame),
    family = if (family == "binomial") binomial("log") else poisson("log"),
    covariance = V, converged = sol$converged, rank = ncol(X),
    fc.diagnostics = sol$diagnostics, fc.bias_model = bias_model,
    information = tt$information, score_rows = tt$score_rows,
    variance_bread = variance_bread, variance_bread_type = variance_bread_type), class = "glm_fc_fit")
  idx <- match("x", names(theta)); se <- sqrt(V[idx, idx]); est <- theta[idx]
  list(est = unname(est), se = unname(se), low = unname(est-qnorm(.975)*se),
    up = unname(est+qnorm(.975)*se), p = if (is.finite(se) && se > 0) 2*pnorm(-abs(est/se)) else NA_real_,
    point.est = theta, vcov = V, vcov.ipd = Vipd, vcov.target = Vtarget,
    fit = fit, converged = sol$converged, fc.diagnostics = sol$diagnostics,
    covariance_type = covariance, variance_status = variance_status,
    variance_bread_type = variance_bread_type,
    variance.method = if (covariance == "model") "fixed_analytic_model_sandwich_first_order" else
      if (account_for_gamma) "joint_empirical_sandwich_first_order" else "fixed_empirical_sandwich_first_order",
    bias_model = bias_model, influence = influence)
}

fit_logbinomial_fc <- function(data, w = NULL, ...) fit_logglm_fc(data, w,
  family = "binomial", bias_model = "bernoulli", covariance = "model", ...)
fit_robust_lb_fc <- function(data, w = NULL, ...) fit_logglm_fc(data, w,
  family = "binomial", bias_model = "bernoulli", covariance = "sandwich", ...)
fit_lp_fc <- function(data, w = NULL, bias_model = "poisson", ...) fit_logglm_fc(data, w,
  family = "poisson", bias_model = bias_model, covariance = "model", ...)
fit_robust_lp_fc <- function(data, w = NULL, bias_model = "poisson", ...) fit_logglm_fc(data, w,
  family = "poisson", bias_model = bias_model, covariance = "sandwich", ...)

fit_anchored_glm_fc <- function(data, w, model_covariates, effect_modifiers,
                               fit_control, use_maic = FALSE, joint_args = list()) {
  base <- list(data = data, w = w, model_covariates = model_covariates,
    effect_modifiers = effect_modifiers, maxit = fit_control$fc_max_step,
    tolerance = fit_control$fc_thres)
  specs <- list(c("binomial", "bernoulli"), c("poisson", "poisson"), c("poisson", "bernoulli"))
  out <- setNames(vector("list", 6L), anchored_glm_fc_method_names())
  for (j in seq_along(specs)) {
    ordinary <- safe_try(do.call(fit_logglm_fc, c(base, list(family = specs[[j]][1],
      bias_model = specs[[j]][2], covariance = "model"))), names(out)[2*j-1L])
    robust <- ordinary
    if (!is.null(ordinary) && isTRUE(ordinary$converged)) {
      robust <- safe_try(do.call(fit_logglm_fc, c(base, list(family = specs[[j]][1],
        bias_model = specs[[j]][2], covariance = "sandwich", start = ordinary$point.est,
        account_for_gamma = use_maic, account_for_target_mean = use_maic), joint_args)), names(out)[2*j])
    } else if (!is.null(robust)) {
      robust$covariance_type <- "sandwich"
      robust$variance.method <- if (use_maic) "joint_empirical_sandwich_first_order" else
        "fixed_empirical_sandwich_first_order"
    }
    out[2*j-1L] <- list(ordinary); out[2*j] <- list(robust)
  }
  out
}
