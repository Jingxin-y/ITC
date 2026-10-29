## BRM、log-binomial、log-Poisson 和联合 sandwich 回归估计。
## 原代码：compare_RR_anchored.R

fit_logbinomial <- function(data, w = NULL, label = "log-binomial",
                            model_covariates = NULL, effect_modifiers = NULL, maxit = 200L) {
  dat <- add_glm_weights(data, w)
  w <- dat$.glm_weight

  tryCatch(
    glm(
      weighted_glm_formula(dat, model_covariates, effect_modifiers),
      family = binomial(link = "log"),
      control = glm.control(maxit = maxit),
      data = dat,
      weights = .glm_weight,
      start = log_link_start(data, w, max_prob = 0.8,
                             model_covariates = model_covariates, effect_modifiers = effect_modifiers)
    ),
    error = function(e) {
      warning(label, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

fit_lp <- function(data, w = NULL, model_covariates = NULL, effect_modifiers = NULL, maxit = 200L) {
  dat <- add_glm_weights(data, w)
  w <- dat$.glm_weight

  glm(
    weighted_glm_formula(dat, model_covariates, effect_modifiers),
    family = poisson(link = "log"),
    control = glm.control(maxit = maxit),
    data = dat,
    weights = .glm_weight,
    start = log_link_start(data, w, max_prob = 0.95,
                           model_covariates = model_covariates, effect_modifiers = effect_modifiers)
  )
}

fit_robust_lp <- function(data, w = NULL, hc_type = "HC0",
                          model_covariates = NULL, effect_modifiers = NULL,
                          return_full = length(effect_modifiers) > 0L, maxit = 200L) {
  if (is.null(w)) w <- rep(1, nrow(data))

  fit <- fit_lp(data, w, model_covariates, effect_modifiers, maxit = maxit)

  est <- unname(coef(fit)["x"])
  covariance <- sandwich::vcovHC(fit, type = hc_type)
  se <- sqrt(unname(covariance["x", "x"]))

  out <- list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = 2 * (1 - pnorm(abs(est / se)))
  )
  if (return_full) {
    out$point.est <- coef(fit)
    out$vcov <- covariance
    out$fit <- fit
  }
  out
}

fit_robust_lp_joint <- function(data, w = NULL, gamma = NULL,
                                target_mean = NULL,
                                target_var = NULL,
                                m_AD = NULL,
                                target_cov = NULL,
                                target_mean_vcov = NULL,
                                account_for_gamma = TRUE,
                                account_for_target_mean = FALSE,
                                weight_covariates = NULL, model_covariates = NULL,
                                effect_modifiers = NULL,
                                return_full = length(effect_modifiers) > 0L, maxit = 200L) {
  ## Joint sandwich for weighted robust log-Poisson. The point estimate is the
  ## same as fit_robust_lp(); the variance additionally propagates uncertainty
  ## from gamma and, optionally, the AD target mean.
  if (is.null(w)) w <- rep(1, nrow(data))
  w <- as.numeric(w)
  if (length(w) != nrow(data)) {
    stop("Length of w must equal nrow(data).")
  }
  if (isTRUE(account_for_target_mean) && !isTRUE(account_for_gamma)) {
    stop("To account for target_mean uncertainty, account_for_gamma should be TRUE.")
  }

  fit <- fit_lp(data, w, model_covariates, effect_modifiers, maxit = maxit)
  beta.hat <- stats::coef(fit)
  coef.idx <- match("x", names(beta.hat))
  est <- unname(beta.hat[coef.idx])

  X <- stats::model.matrix(fit)
  n <- nrow(data)
  if (nrow(X) != n) {
    stop("RLP joint sandwich requires glm to keep all rows.")
  }

  y <- as.numeric(data$y)
  mu <- as.numeric(stats::fitted(fit))
  resid <- y - mu

  psi.beta <- sweep(X, 1, w * resid, "*")
  A.beta <- -crossprod(X, sweep(X, 1, w * mu, "*")) / n
  A.beta.inv <- MASS::ginv(A.beta)

  IF.beta <- -psi.beta %*% t(A.beta.inv)
  var.mu <- matrix(0, nrow = ncol(X), ncol = ncol(X))

  weight_covariates <- weight_covariate_names(data, weight_covariates, gamma)
  Z <- covariate_matrix(data, weight_covariates)
  q <- ncol(Z)

  if (isTRUE(account_for_gamma) && q > 0) {
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
    if (!is.null(gamma)) {
      if (!is.null(names(gamma))) gamma <- gamma[weight_covariates]
      gamma <- as.numeric(gamma)
      if (length(gamma) == 1 && q > 1) gamma <- rep(gamma, q)
    }
    if (!is.null(gamma) && length(gamma) != q) {
      stop("gamma length must match the number of covariates.")
    }

    D <- sweep(Z, 2, target_mean, "-")
    Dbar.w <- colSums(w * D) / sum(w)
    Hgamma <- sweep(D, 2, Dbar.w, "-")

    A.gamma <- crossprod(D, w * Hgamma) / n
    A.gamma.inv <- MASS::ginv(A.gamma)

    A.beta.gamma <- crossprod(
      X,
      sweep(Hgamma, 1, w * resid, "*")
    ) / n

    psi.gamma <- sweep(D, 1, w, "*")
    IF.gamma <- -psi.gamma %*% t(A.gamma.inv)

    IF.beta <- -(
      psi.beta + IF.gamma %*% t(A.beta.gamma)
    ) %*% t(A.beta.inv)

    if (isTRUE(account_for_target_mean)) {
      V.mu <- make_target_mean_vcov(
        target_var = target_var,
        m_AD = m_AD,
        target_cov = target_cov,
        target_mean_vcov = target_mean_vcov
      )

      if (!all(dim(V.mu) == c(q, q))) {
        stop("Dimension mismatch: Var(target_mean_hat) must be q by q.")
      }

      dgamma.dmu <- A.gamma.inv %*% (mean(w) * diag(q))
      dbeta.dmu <- -A.beta.inv %*% A.beta.gamma %*% dgamma.dmu
      var.mu <- dbeta.dmu %*% V.mu %*% t(dbeta.dmu)
    }
  }

  vcov.ipd <- crossprod(IF.beta) / n^2
  vcov.total <- vcov.ipd + var.mu
  se <- sqrt(as.numeric(vcov.total[coef.idx, coef.idx]))

  out <- list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = if (is.finite(se) && se > 0) {
      2 * (1 - pnorm(abs(est / se)))
    } else {
      NA_real_
    },
    var.ipd = as.numeric(vcov.ipd[coef.idx, coef.idx]),
    var.mu = as.numeric(var.mu[coef.idx, coef.idx]),
    vcov = vcov.total
  )
  if (return_full) {
    dimnames(out$vcov) <- list(names(beta.hat), names(beta.hat))
    out$point.est <- beta.hat
    out$fit <- fit
    out$vcov.ipd <- vcov.ipd
    out$vcov.target <- var.mu
  }
  out
}

## 固定权重的稳健 log-binomial：与 fit_robust_lp() 一样直接调用 vcovHC()。
fit_robust_lb <- function(data, w = NULL, model_covariates = NULL,
                          effect_modifiers = NULL,
                          return_full = length(effect_modifiers) > 0L,
                          hc_type = "HC0", maxit = 200L) {
  if (is.null(w)) w <- rep(1, nrow(data))
  w <- as.numeric(w)

  fit <- fit_logbinomial(data, w, model_covariates = model_covariates,
                         effect_modifiers = effect_modifiers, maxit = maxit)
  if (is.null(fit)) return(NULL)

  beta.hat <- stats::coef(fit)
  coef.idx <- match("x", names(beta.hat))
  est <- unname(beta.hat[coef.idx])

  mu <- as.numeric(stats::fitted(fit))
  if (!isTRUE(fit$converged) || any(!is.finite(beta.hat)) || any(mu <= 0 | mu >= 1)) {
    stop("RLB sandwich requires a converged interior log-binomial fit.")
  }
  covariance <- sandwich::vcovHC(fit, type = hc_type)
  se <- sqrt(as.numeric(covariance[coef.idx, coef.idx]))

  out <- list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = if (is.finite(se) && se > 0) {
      2 * (1 - pnorm(abs(est / se)))
    } else {
      NA_real_
    },
    var.ipd = as.numeric(covariance[coef.idx, coef.idx]),
    var.mu = 0,
    vcov = covariance
  )
  if (return_full) {
    dimnames(out$vcov) <- list(names(beta.hat), names(beta.hat))
    out$point.est <- beta.hat
    out$fit <- fit
    out$vcov.ipd <- covariance
    out$vcov.target <- matrix(0, length(beta.hat), length(beta.hat))
  }
  out
}

fit_robust_lb_joint <- function(data, w = NULL, gamma = NULL,
                                target_mean = NULL,
                                target_var = NULL,
                                m_AD = NULL,
                                target_cov = NULL,
                                target_mean_vcov = NULL,
                                account_for_gamma = TRUE,
                                account_for_target_mean = FALSE,
                                weight_covariates = NULL, model_covariates = NULL,
                                effect_modifiers = NULL,
                                return_full = length(effect_modifiers) > 0L, maxit = 200L) {
  ## 与 fit_robust_lp_joint() 保持相同的计算顺序：
  ## 拟合 -> beta 分数/导数 -> gamma 联合方差 -> AD 均值误差 -> 返回结果。
  ## 点估计与 fit_robust_lb() 相同；方差额外传播 gamma 及可选的 AD 均值误差。
  ## 本联合版仍使用观测导数；固定权重版 vcovHC() 使用 GLM 的工作 Fisher bread。
  if (is.null(w)) w <- rep(1, nrow(data))
  w <- as.numeric(w)
  if (length(w) != nrow(data)) {
    stop("Length of w must equal nrow(data).")
  }
  if (isTRUE(account_for_target_mean) && !isTRUE(account_for_gamma)) {
    stop("To account for target_mean uncertainty, account_for_gamma should be TRUE.")
  }

  fit <- fit_logbinomial(data, w, model_covariates = model_covariates,
                         effect_modifiers = effect_modifiers, maxit = maxit)
  if (is.null(fit)) return(NULL)
  beta.hat <- stats::coef(fit)
  coef.idx <- match("x", names(beta.hat))
  est <- unname(beta.hat[coef.idx])

  X <- stats::model.matrix(fit)
  n <- nrow(data)
  if (nrow(X) != n) {
    stop("RLB joint sandwich requires glm to keep all rows.")
  }

  y <- as.numeric(data$y)
  mu <- as.numeric(stats::fitted(fit))
  if (!isTRUE(fit$converged) || any(!is.finite(beta.hat)) || any(mu <= 0 | mu >= 1)) {
    stop("RLB sandwich requires a converged interior log-binomial fit.")
  }
  ## LB 分数和观测导数；相应 LP 位置使用 resid=y-mu 和 w*mu。
  resid <- (y - mu) / (1 - mu)

  psi.beta <- sweep(X, 1, w * resid, "*")
  A.beta <- -crossprod(X, sweep(X, 1, w * mu * (1 - y) / (1 - mu)^2, "*")) / n
  A.beta.inv <- MASS::ginv(A.beta)

  IF.beta <- -psi.beta %*% t(A.beta.inv)
  var.mu <- matrix(0, nrow = ncol(X), ncol = ncol(X))

  weight_covariates <- weight_covariate_names(data, weight_covariates, gamma)
  Z <- covariate_matrix(data, weight_covariates)
  q <- ncol(Z)

  if (isTRUE(account_for_gamma) && q > 0) {
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
    if (!is.null(gamma)) {
      if (!is.null(names(gamma))) gamma <- gamma[weight_covariates]
      gamma <- as.numeric(gamma)
      if (length(gamma) == 1 && q > 1) gamma <- rep(gamma, q)
    }
    if (!is.null(gamma) && length(gamma) != q) {
      stop("gamma length must match the number of covariates.")
    }

    D <- sweep(Z, 2, target_mean, "-")
    Dbar.w <- colSums(w * D) / sum(w)
    Hgamma <- sweep(D, 2, Dbar.w, "-")

    A.gamma <- crossprod(D, w * Hgamma) / n
    A.gamma.inv <- MASS::ginv(A.gamma)

    A.beta.gamma <- crossprod(
      X,
      sweep(Hgamma, 1, w * resid, "*")
    ) / n

    psi.gamma <- sweep(D, 1, w, "*")
    IF.gamma <- -psi.gamma %*% t(A.gamma.inv)

    IF.beta <- -(
      psi.beta + IF.gamma %*% t(A.beta.gamma)
    ) %*% t(A.beta.inv)

    if (isTRUE(account_for_target_mean)) {
      V.mu <- make_target_mean_vcov(
        target_var = target_var,
        m_AD = m_AD,
        target_cov = target_cov,
        target_mean_vcov = target_mean_vcov
      )

      if (!all(dim(V.mu) == c(q, q))) {
        stop("Dimension mismatch: Var(target_mean_hat) must be q by q.")
      }

      dgamma.dmu <- A.gamma.inv %*% (mean(w) * diag(q))
      dbeta.dmu <- -A.beta.inv %*% A.beta.gamma %*% dgamma.dmu
      var.mu <- dbeta.dmu %*% V.mu %*% t(dbeta.dmu)
    }
  }

  vcov.ipd <- crossprod(IF.beta) / n^2
  vcov.total <- vcov.ipd + var.mu
  se <- sqrt(as.numeric(vcov.total[coef.idx, coef.idx]))

  out <- list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = if (is.finite(se) && se > 0) {
      2 * (1 - pnorm(abs(est / se)))
    } else {
      NA_real_
    },
    var.ipd = as.numeric(vcov.ipd[coef.idx, coef.idx]),
    var.mu = as.numeric(var.mu[coef.idx, coef.idx]),
    vcov = vcov.total
  )
  if (return_full) {
    dimnames(out$vcov) <- list(names(beta.hat), names(beta.hat))
    out$point.est <- beta.hat
    out$fit <- fit
    out$vcov.ipd <- vcov.ipd
    out$vcov.target <- var.mu
  }
  out
}
