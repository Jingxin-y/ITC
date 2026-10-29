## 分层 Mantel-Haenszel 估计及其 sandwich 方差。
## 原代码：compare_RR_anchored.R

get.weighted.RR <- function(data, w = NULL, correction = 0.5,
                            strata = NULL, return_details = FALSE,
                            model_covariates = NULL) {
  ## Stratified Mantel-Haenszel estimator for the common log risk ratio.
  ##
  ## For stratum s, let
  ##   a_s = sum w I(x=1,y=1), b_s = sum w I(x=1,y=0),
  ##   c_s = sum w I(x=0,y=1), d_s = sum w I(x=0,y=0),
  ##   n1_s = a_s+b_s, n0_s = c_s+d_s, N_s=n1_s+n0_s.
  ## The MH estimator is
  ##   RR_MH = sum_s a_s n0_s / N_s / sum_s c_s n1_s / N_s.
  ## The SE below is a delta-method / influence-function SE treating w as fixed.

  n <- nrow(data)
  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w)
  if (length(w) != n) stop("Length of w must equal nrow(data).")

  x <- data$x
  y <- data$y

  if (is.null(strata)) {
    Z <- as.matrix(covariate_matrix(data, model_covariates))
    if (ncol(Z) == 0) {
      strata <- factor(rep("all", n))
    } else {
      strata <- do.call(
        interaction,
        c(as.data.frame(Z), list(drop = TRUE, lex.order = TRUE, sep = ":"))
      )
    }
  } else {
    if (length(strata) != n) stop("Length of strata must equal nrow(data).")
    strata <- factor(strata)
  }

  lev <- levels(strata)
  raw <- do.call(rbind, lapply(lev, function(s) {
    idx <- strata == s
    a <- sum(w[idx] * (x[idx] == 1) * (y[idx] == 1))
    b <- sum(w[idx] * (x[idx] == 1) * (y[idx] == 0))
    c <- sum(w[idx] * (x[idx] == 0) * (y[idx] == 1))
    d <- sum(w[idx] * (x[idx] == 0) * (y[idx] == 0))

    data.frame(
      stratum = s,
      a = a,
      b = b,
      c = c,
      d = d,
      n1 = a + b,
      n0 = c + d,
      stringsAsFactors = FALSE
    )
  }))

  ## A stratum contributes to MH only if both treatment arms are present.
  keep <- raw$n1 > 0 & raw$n0 > 0
  raw <- raw[keep, , drop = FALSE]

  make_na <- function() {
    out <- list(est = NA_real_, se = NA_real_,
                low = NA_real_, up = NA_real_, p = NA_real_)
    if (isTRUE(return_details)) {
      out$phi <- rep(NA_real_, n)
      out$grad <- NULL
      out$part <- NULL
    }
    out
  }

  if (nrow(raw) == 0) return(make_na())

  counts <- raw

  mh_num_den <- function(cnt) {
    N <- cnt$n1 + cnt$n0
    c(
      num = sum(cnt$a * cnt$n0 / N),
      den = sum(cnt$c * cnt$n1 / N)
    )
  }

  nd <- mh_num_den(counts)
  used_correction <- FALSE

  ## If the MH numerator or denominator is zero, apply a continuity correction
  ## to the contributing strata. This keeps log(RR_MH) finite in rare-event runs.
  if (!all(is.finite(nd)) || nd["num"] <= 0 || nd["den"] <= 0) {
    counts[, c("a", "b", "c", "d")] <-
      counts[, c("a", "b", "c", "d")] + correction
    counts$n1 <- counts$a + counts$b
    counts$n0 <- counts$c + counts$d
    nd <- mh_num_den(counts)
    used_correction <- TRUE
  }

  if (!all(is.finite(nd)) || nd["num"] <= 0 || nd["den"] <= 0) {
    return(make_na())
  }

  cmh_grad <- function(cnt) {
    N <- cnt$n1 + cnt$n0
    A <- sum(cnt$a * cnt$n0 / N)
    C <- sum(cnt$c * cnt$n1 / N)

    list(
      A = A,
      C = C,
      est = log(A / C),
      ga = (cnt$n0 / N) / A,
      gc = -(cnt$n1 / N) / C,
      gn1 = -(cnt$a * cnt$n0 / N^2) / A -
        (cnt$c * cnt$n0 / N^2) / C,
      gn0 = (cnt$a * cnt$n1 / N^2) / A +
        (cnt$c * cnt$n1 / N^2) / C
    )
  }

  ## Use mean-scale quantities for the influence-function calculation.
  counts.mean <- counts
  counts.mean[, c("a", "b", "c", "d", "n1", "n0")] <-
    counts.mean[, c("a", "b", "c", "d", "n1", "n0")] / n

  raw.mean <- raw
  raw.mean[, c("a", "b", "c", "d", "n1", "n0")] <-
    raw.mean[, c("a", "b", "c", "d", "n1", "n0")] / n

  grad <- cmh_grad(counts.mean)

  kept.strata <- factor(strata, levels = counts$stratum)
  gobs <- rep(0, n)

  for (k in seq_along(counts$stratum)) {
    idx <- !is.na(kept.strata) & kept.strata == counts$stratum[k]
    if (!any(idx)) next

    gobs[idx] <- grad$ga[k]  * (x[idx] == 1) * (y[idx] == 1) +
      grad$gc[k]  * (x[idx] == 0) * (y[idx] == 1) +
      grad$gn1[k] * (x[idx] == 1) +
      grad$gn0[k] * (x[idx] == 0)
  }

  center <- sum(
    grad$ga  * raw.mean$a +
      grad$gc  * raw.mean$c +
      grad$gn1 * raw.mean$n1 +
      grad$gn0 * raw.mean$n0
  )

  phi <- w * gobs - center
  se <- sqrt(sum(phi^2) / n^2)
  est <- grad$est

  out <- list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = if (is.finite(se) && se > 0) {
      2 * (1 - pnorm(abs(est / se)))
    } else {
      NA_real_
    }
  )

  if (isTRUE(return_details)) {
    out$phi <- phi
    out$gobs <- gobs
    out$grad <- grad
    out$part <- list(
      raw = raw,
      counts = counts,
      strata = kept.strata,
      used_correction = used_correction
    )
  }

  out
}

get.weighted.RR.sandwich <- function(data, w, gamma = NULL,
                                     target_mean = NULL,
                                     target_var = NULL,
                                     m_AD = NULL,
                                     target_cov = NULL,
                                     target_mean_vcov = NULL,
                                     account_for_gamma = TRUE,
                                     account_for_target_mean = FALSE,
                                     correction = 0.5,
                                     strata = NULL,
                                     weight_covariates = NULL, model_covariates = NULL) {
  ## Stratified Mantel-Haenszel common log RR with sandwich SE.
  ## The fixed-weight component is the same CMH influence-function SE used in
  ## get.weighted.RR(). If account_for_gamma = TRUE, the SE also accounts for
  ## estimation of the exponential-tilting parameter gamma. If
  ## account_for_target_mean = TRUE, it additionally accounts for uncertainty
  ## in the AD target mean used to estimate gamma.

  fixed <- get.weighted.RR(
    data = data,
    w = w,
    correction = correction,
    strata = strata,
    return_details = TRUE,
    model_covariates = model_covariates
  )

  base_out <- list(
    est = fixed$est,
    se = fixed$se,
    low = fixed$low,
    up = fixed$up,
    p = fixed$p
  )

  if (!is.finite(fixed$est) || is.null(fixed$phi) || is.null(fixed$part)) {
    return(base_out)
  }

  use_gamma <- isTRUE(account_for_gamma)
  use_mu <- isTRUE(account_for_target_mean)

  if (use_mu && !use_gamma) {
    stop("To account for target_mean uncertainty, account_for_gamma should be TRUE.")
  }

  if (!use_gamma) {
    return(base_out)
  }

  n <- nrow(data)
  w <- as.numeric(w)

  weight_covariates <- weight_covariate_names(data, weight_covariates, gamma)
  Z <- covariate_matrix(data, weight_covariates)
  Z <- as.matrix(Z)
  q <- ncol(Z)

  if (q == 0) {
    return(base_out)
  }

  if (is.null(target_mean)) {
    stop("target_mean must be supplied when account_for_gamma = TRUE.")
  }

  all_covariates <- covariate_names(data)
  target_mean <- as.numeric(select_weight_summary(target_mean, weight_covariates, all_covariates))
  target_var <- select_weight_summary(target_var, weight_covariates, all_covariates)
  target_cov <- select_weight_summary(target_cov, weight_covariates, all_covariates)
  target_mean_vcov <- select_weight_summary(target_mean_vcov, weight_covariates, all_covariates)
  if (length(target_mean) == 1 && q > 1) target_mean <- rep(target_mean, q)
  if (length(target_mean) != q) {
    stop("target_mean length must match the number of covariates.")
  }

  ## Exponential-tilting calibration equation:
  ##   mean{ w_gamma(Z_i) (Z_i - target_mean) } = 0.
  ## Because weights are normalized, d w_i / d gamma = w_i H_i, where
  ##   H_i = D_i - weighted.mean(D).
  D <- sweep(Z, 2, target_mean, "-")
  Dbar.w <- colSums(w * D) / sum(w)
  Hgamma <- sweep(D, 2, Dbar.w, "-")

  A.gamma <- crossprod(D, w * Hgamma) / n
  A.gamma.inv <- MASS::ginv(A.gamma)

  ## Derivative of the CMH log RR wrt gamma, holding the empirical
  ## distribution fixed.
  G.gamma <- colSums((w * fixed$gobs) * Hgamma) / n

  ## IF for gamma from the IPD calibration equation.
  psi.gamma <- w * D
  IF.gamma <- -psi.gamma %*% t(A.gamma.inv)

  phi.ipd <- fixed$phi + as.vector(IF.gamma %*% G.gamma)
  var.ipd <- sum(phi.ipd^2) / n^2

  var.mu <- 0

  if (use_mu) {
    V.mu <- make_target_mean_vcov(
      target_var = target_var,
      m_AD = m_AD,
      target_cov = target_cov,
      target_mean_vcov = target_mean_vcov
    )

    if (!all(dim(V.mu) == c(q, q))) {
      stop("Dimension mismatch: Var(target_mean_hat) must be q by q.")
    }

    ## Since normalized weights are invariant to a common shift in
    ## Z - target_mean, target_mean affects the estimator through gamma.
    ## From U_gamma(gamma, mu)=0,
    ##   d gamma / d mu = A.gamma^{-1} mean(w) I.
    dgamma.dmu <- A.gamma.inv %*% (mean(w) * diag(q))
    G.mu <- as.vector(t(G.gamma) %*% dgamma.dmu)
    var.mu <- as.numeric(t(G.mu) %*% V.mu %*% G.mu)
  }

  var.total <- var.ipd + var.mu
  se <- sqrt(var.total)
  est <- fixed$est

  list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = if (is.finite(se) && se > 0) {
      2 * (1 - pnorm(abs(est / se)))
    } else {
      NA_real_
    },
    var.ipd = var.ipd,
    var.mu = var.mu
  )
}
