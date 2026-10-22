rd_validate_seeds <- function(seeds) {
  if (!is.numeric(seeds) || !length(seeds) || any(!is.finite(seeds)) ||
      any(seeds < 0 | seeds > .Machine$integer.max | seeds != floor(seeds)) || anyDuplicated(seeds))
    stop("seeds must be distinct nonnegative R integers.")
}
rd_validate_sizes <- function(n_AC, n_BC) {
  for (n in list(n_AC, n_BC)) if (!is.numeric(n) || length(n) != 1L ||
      !is.finite(n) || n < 2 || n > .Machine$integer.max || n %% 2 != 0)
    stop("n_AC/n_BC must be positive even total trial sizes.")
}
rd_weights <- function(w, n) {
  if (is.null(w)) w <- rep(1, n)
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w)) || any(w < 0) || !any(w > 0))
    stop("Weights must be finite, nonnegative, and have positive total.")
  w <- w / max(w); w / mean(w)
}
rd_method_plan <- function() data.frame(
  method = c("brm", "brm-BC", "brm-FC", "brm-FC-BC", "IB", "RIB", "RLPM", "MN_eff", "RD_MAIC_sandwich"),
  fit_method = c("brm", "brm", "brm-FC", "brm-FC", "IB", "IB", "RLPM", "MN_eff", "RD_MAIC_sandwich"),
  ci_type = c("Wald", "BC", "Wald", "BC", "Wald", "Wald", "Wald", "MN-ESS", "Wald"),
  stringsAsFactors = FALSE)

## Numerically stable RD/OP map, including small risks and zero RD/OP.
## This is for estimation/CI only; the DGM never calls it.
rd_probabilities <- function(theta, phi) {
  if (!is.numeric(theta) || !is.numeric(phi) || is.complex(theta) || is.complex(phi) ||
      any(!is.finite(theta)) || any(!is.finite(phi)))
    stop("RD/OP predictors must be finite real numeric values.")
  theta <- as.numeric(theta); phi <- as.numeric(phi)
  if (!length(theta) && !length(phi)) return(list(p0 = numeric(), p1 = numeric()))
  n <- max(length(theta), length(phi))
  if (!length(theta) || !length(phi) || !all(c(length(theta), length(phi)) %in% c(1L, n)))
    stop("RD/OP predictors must have equal lengths or be scalar.")
  theta <- rep_len(theta, n); phi <- rep_len(phi, n)
  d <- tanh(theta); t <- abs(d); hi <- phi >= 0; m <- numeric(n)
  ## Rationalized quadratic discriminant: hypot(t*(1-u), 2*sqrt(u)).
  ## Rescaling avoids squaring tiny values; never exponentiate positive phi.
  u <- exp(-abs(phi)); a <- t * (-expm1(-abs(phi))); b <- 2 * exp(-abs(phi) / 2)
  scale <- pmax(a, b); root <- numeric(n); active <- scale > 0
  root[active] <- scale[active] * sqrt((a[active] / scale[active])^2 + (b[active] / scale[active])^2)
  zero <- t == 0
  m[zero] <- plogis(phi[zero] / 2)
  pos <- hi & !zero
  m[pos] <- 2 * (1 - t[pos]) / (2 - t[pos] + t[pos] * u[pos] + root[pos])
  neg <- !hi & !zero & t < 1
  ## u may underflow although u/t is representable: evaluate the numerator in log space.
  den <- u[neg] * (2 - t[neg]) + t[neg] + root[neg]
  m[neg] <- exp(log(2) + phi[neg] + log1p(-t[neg]) - log(den))
  ## When tanh(theta) rounds to +/-1 the compatible risks round to 0/1.
  ## This is machine precision, not the legacy +/-12 OP threshold approximation.
  list(p0 = m + ifelse(d < 0, t, 0), p1 = m + ifelse(d > 0, t, 0))
}
marginal_rd_target <- function(pars, va, vb, weight = NULL, gradient = FALSE) {
  va <- as.matrix(va); vb <- as.matrix(vb); pa <- ncol(va)
  if (!is.numeric(pars) || length(pars) != pa + ncol(vb) || nrow(va) != nrow(vb) ||
      !nrow(va) || pa < 1 || ncol(vb) < 1 || any(!is.finite(c(pars, va, vb))))
    stop("Invalid RD target parameters/designs.")
  w <- rd_weights(weight, nrow(va)); omega <- w / sum(w)
  theta <- drop(va %*% pars[seq_len(pa)]); d <- tanh(theta)
  value <- sum(omega * d)
  if (!gradient) return(value)
  p <- rd_probabilities(theta, drop(vb %*% pars[pa + seq_len(ncol(vb))]))
  list(value = value, gradient = c(colSums(va * (omega * (1 - d^2))), rep(0, ncol(vb))),
    risk.C = sum(omega * p$p0), risk.A = sum(omega * p$p1))
}

rd_glm_design <- function(data, model_covariates = c("h_E", "h_P"), effect_modifiers = "h_E") {
  f <- weighted_glm_formula(data, model_covariates, effect_modifiers)
  X <- model.matrix(f, data); d0 <- d1 <- data; d0$x <- 0; d1$x <- 1
  list(formula = f, X = X, X0 = model.matrix(f, d0), X1 = model.matrix(f, d1))
}
rd_fit_identity <- function(data, w, design, maxit = 200L) {
  X <- design$X; w <- rd_weights(w, nrow(X))
  if (qr(X[w > 0, , drop = FALSE])$rank != ncol(X)) stop("Identity design is not identifiable.")
  start <- numeric(ncol(X)); start[1] <- weighted.mean(data$y, w)
  start[1] <- min(.99, max(.01, start[1]))
  fit <- glm.fit(X, data$y, weights = w, family = binomial("identity"), start = start,
    control = glm.control(maxit = maxit, epsilon = 1e-9))
  if (!isTRUE(fit$converged)) stop("Identity-binomial did not converge.")
  p <- fit$fitted.values
  if (any(p <= 0 | p >= 1) || any(!is.finite(fit$coefficients))) stop("Identity fit is not interior.")
  v <- p * (1 - p); score <- X * ((data$y - p) / v)
  H <- crossprod(X, X * (w * (data$y / p^2 + (1 - data$y) / (1 - p)^2))) / nrow(X)
  info <- crossprod(X, X * (w / v)); Ii <- solve(info)
  covariance <- Ii %*% crossprod(X, X * (w^2 / v)) %*% Ii
  robust <- solve(H) %*% crossprod(score * w) %*% t(solve(H)) / nrow(X)^2
  list(coefficients = fit$coefficients, converged = TRUE, fit = fit, design = design,
    covariance = covariance, robust_covariance = robust, score = score, H = H,
    influence_score = score, factual_risks_valid = TRUE)
}
rd_fit_lpm <- function(data, w, design) {
  X <- design$X; w <- rd_weights(w, nrow(X)); n <- nrow(X)
  if (qr(X[w > 0, , drop = FALSE])$rank != ncol(X)) stop("LPM design is not identifiable.")
  fit <- lm.wfit(X, data$y, w); e <- fit$residuals; Ai <- solve(crossprod(X, X * w))
  h <- w * rowSums((X %*% Ai) * X)
  if (any(1 - h <= 1e-10)) stop("HC3 undefined at unit leverage.")
  hc_score <- X * (e / (1 - h))
  V <- Ai %*% crossprod(hc_score * w) %*% Ai
  list(coefficients = fit$coefficients, converged = TRUE, fit = fit, design = design,
    covariance = V, robust_covariance = V, score = X * e,
    influence_score = hc_score, H = crossprod(X, X * w) / n, leverage = h,
    factual_risks_valid = all(fit$fitted.values >= 0 & fit$fitted.values <= 1))
}
rd_regression_terms <- function(fit, robust = FALSE) {
  X0 <- fit$design$X0; X1 <- fit$design$X1
  list(p0 = drop(X0 %*% fit$coefficients), p1 = drop(X1 %*% fit$coefficients),
    dp0 = X0, dp1 = X1, score = fit$score, influence_score = fit$influence_score,
    H = fit$H, covariance = if (robust) fit$robust_covariance else fit$covariance)
}
rd_brm_terms <- function(fit, y, x, w) {
  va <- fit$va; vb <- fit$vb; pa <- ncol(va); n <- length(y)
  a <- fit$point.est[seq_len(pa)]; b <- fit$point.est[pa + seq_len(ncol(vb))]
  tt <- rd_solver_derivatives(a, b, x, va, vb)
  D <- cbind(va * tt$d.a, vb * tt$d.b); p <- tt$p; v <- p * (1 - p)
  residual <- (y - p) / v
  hp <- rbind(cbind(crossprod(va, va * (w * residual * tt$d.aa)),
                   crossprod(va, vb * (w * residual * tt$d.ab))),
              cbind(crossprod(vb, va * (w * residual * tt$d.ab)),
                    crossprod(vb, vb * (w * residual * tt$d.bb))))
  H <- (crossprod(D, D * (w * (y / p^2 + (1 - y) / (1 - p)^2))) - hp) / n
  z <- rd_probabilities(drop(va %*% a), drop(vb %*% b))
  s0 <- z$p0 * (1 - z$p0); s1 <- z$p1 * (1 - z$p1); denom <- s0 + s1
  da0 <- -s0 / denom * (1 - tanh(drop(va %*% a))^2)
  da1 <- s1 / denom * (1 - tanh(drop(va %*% a))^2)
  db <- s0 * s1 / denom
  list(p0 = z$p0, p1 = z$p1, dp0 = cbind(va * da0, vb * db),
    dp1 = cbind(va * da1, vb * db), score = D * residual,
    influence_score = D * residual, H = H, covariance = fit$cov)
}

## Joint MAIC/standardization/AD-mean influence calculation on the RD scale.
## FC uses the same first-order unadjusted-score convention as existing RR.
rd_delta <- function(terms, v, wfit, Z = NULL, mu = NULL, Vmu = NULL, joint = FALSE) {
  n <- length(v); v <- rd_weights(v, n); wfit <- rd_weights(wfit, n)
  d <- terms$p1 - terms$p0; estimate <- mean(v * d)
  gradient <- colMeans(v * (terms$dp1 - terms$dp0))
  empirical <- v * (d - estimate); variance_target <- 0
  if (!joint) variance_ipd <- drop(crossprod(gradient, terms$covariance %*% gradient)) else {
    D <- sweep(Z, 2, mu, "-"); Hg <- sweep(D, 2, colMeans(v * D), "-")
    Ag <- crossprod(D, v * Hg) / n
    Ai <- if (ncol(Ag)) MASS::ginv(Ag) else matrix(numeric(), 0, 0)
    Hi <- MASS::ginv(terms$H)
    IFg <- -(v * D) %*% t(Ai)
    G <- crossprod(terms$score, wfit * Hg) / n
    IFtheta <- (wfit * terms$influence_score + IFg %*% t(G)) %*% t(Hi)
    gg <- colMeans(Hg * empirical)
    influence <- drop(empirical + IFtheta %*% gradient + IFg %*% gg)
    tg <- drop(t(Hi %*% G %*% Ai) %*% gradient + t(Ai) %*% gg)
    variance_ipd <- sum(influence^2) / n^2
    if (length(tg)) variance_target <- drop(crossprod(tg, Vmu %*% tg))
  }
  out <- rd_wald(estimate, variance_ipd + variance_target)
  c(out, list(risk_C = mean(v * terms$p0), risk_A = mean(v * terms$p1),
    variance_ipd = variance_ipd, variance_target = variance_target,
    counterfactual_risks_valid = all(terms$p0 >= 0 & terms$p0 <= 1 & terms$p1 >= 0 & terms$p1 <= 1)))
}
rd_wald <- function(estimate, variance, level = .95) {
  se <- if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
  z <- qnorm((1 + level) / 2)
  p <- if (is.finite(se) && se > 0) 2 * pnorm(-abs(estimate / se)) else
    if (is.finite(se) && se == 0) as.numeric(estimate == 0) else NA_real_
  list(estimate = estimate, se = se, lower = estimate - z * se, upper = estimate + z * se, p.value = p)
}

rd_mn_score <- function(p1, n1, p0, n0, difference) {
  if (abs(difference) > 1 || min(n1, n0) <= 0 || n1 + n0 <= 1) stop("Invalid MN inputs.")
  lower <- max(0, -difference); upper <- min(1, 1 - difference)
  loglik <- function(q0) {
    qs <- c(q0 + difference, q0); successes <- c(n1 * p1, n0 * p0)
    failures <- c(n1 * (1 - p1), n0 * (1 - p0))
    a <- ifelse(successes == 0, 0, successes * log(qs))
    b <- ifelse(failures == 0, 0, failures * log1p(-qs))
    sum(a + b)
  }
  if (upper - lower < 1e-14) q0 <- lower else {
    mid <- optimize(function(q) -loglik(q), c(lower, upper), tol = 1e-12)$minimum
    candidates <- c(lower, mid, upper); q0 <- candidates[which.max(vapply(candidates, loglik, numeric(1)))]
  }
  variance <- ((q0 + difference) * (1 - q0 - difference) / n1 + q0 * (1 - q0) / n0) *
    (n1 + n0) / (n1 + n0 - 1)
  numerator <- (p1 - p0 - difference)^2
  if (variance <= 1e-30) return(if (numerator <= 1e-24) 0 else Inf)
  numerator / variance
}
rd_mn_interval <- function(p1, n1, p0, n0, level = .95) {
  if (any(!is.finite(c(p1, p0, n1, n0, level))) || min(p1, p0) < 0 || max(p1, p0) > 1 ||
      min(n1, n0) <= 0 || n1 + n0 <= 1 || level <= 0 || level >= 1) stop("Invalid MN inputs.")
  estimate <- p1 - p0; cutoff <- qchisq(level, 1)
  root <- function(boundary) {
    if (abs(boundary - estimate) < 1e-12) return(boundary)
    if (rd_mn_score(p1, n1, p0, n0, boundary) <= cutoff) return(boundary)
    lo <- min(boundary, estimate); hi <- max(boundary, estimate)
    uniroot(function(d) rd_mn_score(p1, n1, p0, n0, d) - cutoff, c(lo, hi), tol = 1e-9)$root
  }
  low <- root(-1); up <- root(1)
  list(estimate = estimate, lower = low, upper = up,
    p.value = pchisq(rd_mn_score(p1, n1, p0, n0, 0), 1, lower.tail = FALSE),
    se = NA_real_, se_from_ci_width = (up - low) / (2 * qnorm((1 + level) / 2)),
    risk_C = p0, risk_A = p1, n_eff_C = n0, n_eff_A = n1,
    approximation = "Kish effective sample size MN; fixed MAIC weights; effective counts may be fractional")
}
rd_weighted_arms <- function(y, x, w) {
  w <- rd_weights(w, length(y)); ps <- ess <- numeric(2)
  for (a in 0:1) {
    idx <- x == a
    if (!any(idx) || sum(w[idx]) <= 0) stop("Both arms need positive weight.")
    ps[a + 1] <- weighted.mean(y[idx], w[idx]); ess[a + 1] <- sum(w[idx])^2 / sum(w[idx]^2)
  }
  list(p0 = ps[1], p1 = ps[2], n0 = ess[1], n1 = ess[2], w = w)
}
rd_weighted_sandwich <- function(y, x, w, Z = NULL, mu = NULL, Vmu = NULL, joint = FALSE) {
  z <- rd_weighted_arms(y, x, w); w <- z$w; n <- length(y)
  phi <- w * (x * (y - z$p1) / mean(w * x) - (1 - x) * (y - z$p0) / mean(w * (1 - x)))
  vt <- 0
  if (joint && ncol(Z)) {
    D <- sweep(Z, 2, mu, "-"); Hg <- sweep(D, 2, colMeans(w * D), "-")
    Ai <- MASS::ginv(crossprod(D, w * Hg) / n)
    dg <- colMeans(Hg * phi); IFg <- -(w * D) %*% t(Ai)
    phi <- phi + drop(IFg %*% dg)
    tg <- drop(t(Ai) %*% dg); vt <- drop(crossprod(tg, Vmu %*% tg))
  }
  vi <- sum(phi^2) / n^2
  c(rd_wald(z$p1 - z$p0, vi + vt), list(risk_C = z$p0, risk_A = z$p1,
    variance_ipd = vi, variance_target = vt, n_eff_C = z$n0, n_eff_A = z$n1))
}
