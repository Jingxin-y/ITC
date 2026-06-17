## Standalone anchored-only helpers for the ITC RR simulation.
##
## Required before sourcing this file:
##   getProbScalarRR.R, MyFunc.R, 1.1_MLE_Point.R, 1.2_MLE_Var.R,
##   1.2_MLE_Var_joint.R, 1_CallMLE.R, CI_exact_fast_small.R,
##   and RcppExports.R when exact() uses mat_vec_mul().

safe_try <- function(expr, label) {
  tryCatch(
    expr,
    error = function(e) {
      warning(label, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

first_value <- function(x) {
  if (length(x) == 0) NA_real_ else x[1]
}

optional_arg <- function(x) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) {
    return(NULL)
  }
  x
}

make_na_est <- function() {
  list(point.est = NA_real_, se.est = NA_real_,
       conf.lower = NA_real_, conf.upper = NA_real_, p.value = NA_real_)
}

make_na_exact <- function() {
  list(low = NA_real_, up = NA_real_, p = NA_real_)
}

stat_value <- function(x, field) {
  if (is.null(x) || is.null(x[[field]])) return(NA_real_)
  first_value(x[[field]])
}



get.aggre <- function(data) {
  x <- data$data$x
  Z <- covariate_matrix(data$data)
  Na0 <- data$count[1]
  Na1 <- data$count[2]
  m.AD <- nrow(Z)
  q <- ncol(Z)

  mean0 <- colMeans(Z[x == 0, , drop = FALSE])
  mean1 <- colMeans(Z[x == 1, , drop = FALSE])
  anchored.mean <- (Na0 * mean0 + Na1 * mean1) / (Na0 + Na1)
  anchored.cov <- if (m.AD > 1) {
    as.matrix(stats::cov(Z))
  } else {
    matrix(NA_real_, nrow = q, ncol = q,
           dimnames = list(colnames(Z), colnames(Z)))
  }
  anchored.var <- diag(anchored.cov)

  list(
    anchored.mean.x = anchored.mean,
    anchored.var.x = anchored.var,
    anchored.cov.x = anchored.cov,
    anchored.n = m.AD
  )
}

get.weight <- function(data, target_mean) {
  Z <- covariate_matrix(data$data)
  q <- ncol(Z)

  if (q == 0) {
    return(list(weight = rep(1, nrow(data$data)), ga = numeric(0)))
  }

  target_mean <- as.numeric(target_mean)
  if (length(target_mean) == 1 && q > 1) target_mean <- rep(target_mean, q)
  if (length(target_mean) != q) {
    stop("target_mean length must match the number of covariates.")
  }

  X <- sweep(Z, 2, target_mean, "-")
  Q <- function(alpha, X) log_sum_exp(X %*% alpha)
  Q.grad <- function(alpha, X) {
    eta <- as.vector(X %*% alpha)
    w <- exp(eta - max(eta))
    as.vector(crossprod(X, w / sum(w)))
  }

  fit <- optim(par = rep(0, q), fn = Q, gr = Q.grad, X = X,
               method = "BFGS", control = list(maxit = 1000))
  if (fit$convergence != 0) {
    warning("Weight optimization did not fully converge: ",
            fit$message, call. = FALSE)
  }

  w <- normalize_exp_weights(X %*% fit$par)
  list(weight = w, ga = fit$par)
}

get.estimate <- function(group, data, weight, gamma.ipd = NULL,
                         target_mean = 0,
                         target_var = NULL,
                         m_AD = NULL,
                         target_cov = NULL,
                         target_mean_vcov = NULL) {
  y <- data$data$y
  x <- data$data$x
  va <- as.matrix(data$data$v.1, ncol = 1)
  vb <- cbind(data$data$v.1, covariate_matrix(data$data))

  pa <- 1
  pb <- ncol(vb)
  if (is.null(gamma.ipd)) gamma.ipd <- rep(0, pb - 1)
  if (length(gamma.ipd) == 1 && pb > 2) gamma.ipd <- rep(gamma.ipd, pb - 1)
  if (length(target_mean) == 1 && pb > 2) target_mean <- rep(target_mean, pb - 1)

  MLEst(
    group, "RR", y, x, va, vb, weight, gamma.ipd,
    max.step = min(pa * 20, 1000),
    thres = 1e-6,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb),
    pa = pa,
    pb = pb,
    target_mean = target_mean,
    target_var = target_var,
    m_AD = m_AD,
    target_cov = target_cov,
    target_mean_vcov = target_mean_vcov
  )
}

log_link_start <- function(data, w = NULL, correction = 0.5,
                           max_prob = 0.8) {
  if (is.null(w)) w <- rep(1, nrow(data))
  terms <- analysis_terms(data)
  start <- rep(0, length(terms))
  names(start) <- terms

  x <- data$x
  y <- data$y
  idx1 <- x == 1
  idx0 <- x == 0

  p1 <- (sum(w[idx1] * y[idx1]) + correction) /
    (sum(w[idx1]) + 2 * correction)
  p0 <- (sum(w[idx0] * y[idx0]) + correction) /
    (sum(w[idx0]) + 2 * correction)

  p1 <- pmin(pmax(p1, 1e-8), max_prob)
  p0 <- pmin(pmax(p0, 1e-8), max_prob)

  start["x"] <- log(p1 / p0)
  start["v.1"] <- log(p0)
  start[is.na(start)] <- 0
  start
}

fit_logbinomial <- function(data, w = NULL, label = "log-binomial") {
  if (is.null(w)) w <- rep(1, nrow(data))

  tryCatch(
    glm(
      weighted_glm_formula(data),
      family = binomial(link = "log"),
      data = data,
      weights = w,
      start = log_link_start(data, w, max_prob = 0.8)
    ),
    error = function(e) {
      warning(label, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

get_glm_stat <- function(fit, stat) {
  if (is.null(fit)) return(NA_real_)

  val <- tryCatch({
    switch(
      stat,
      estimate = unname(stats::coef(fit)[1]),
      se = unname(summary(fit)$coefficients[1, 2]),
      low = unname(confint.default(fit, level = 0.95)[1, 1]),
      up = unname(confint.default(fit, level = 0.95)[1, 2]),
      p = unname(summary(fit)$coefficients[1, 4]),
      NA_real_
    )
  }, error = function(e) NA_real_)

  first_value(val)
}

fit_lp <- function(data, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(data))

  glm(
    weighted_glm_formula(data),
    family = poisson(link = "log"),
    data = data,
    weights = w,
    start = log_link_start(data, w, max_prob = 0.95)
  )
}

fit_robust_lp <- function(data, w = NULL, hc_type = "HC0") {
  if (is.null(w)) w <- rep(1, nrow(data))

  fit <- fit_lp(data, w)

  est <- unname(coef(fit)[1])
  se <- sqrt(sandwich::vcovHC(fit, type = hc_type)[1, 1])

  list(
    est = est,
    se = se,
    low = est - 1.96 * se,
    up = est + 1.96 * se,
    p = 2 * (1 - pnorm(abs(est / se)))
  )
}

get.weighted.RR <- function(data, w, correction = 0.5) {
  x <- data$x
  y <- data$y

  w1 <- w[x == 1]
  y1 <- y[x == 1]
  w0 <- w[x == 0]
  y0 <- y[x == 0]

  p1 <- sum(w1 * y1) / sum(w1)
  p0 <- sum(w0 * y0) / sum(w0)

  if (p1 <= 0 || p0 <= 0 || p1 >= 1 || p0 >= 1) {
    p1 <- (sum(w1 * y1) + correction) / (sum(w1) + 2 * correction)
    p0 <- (sum(w0 * y0) + correction) / (sum(w0) + 2 * correction)
  }

  logrr <- log(p1 / p0)
  n1.eff <- sum(w1)^2 / sum(w1^2)
  n0.eff <- sum(w0)^2 / sum(w0^2)
  se <- sqrt((1 - p1) / (n1.eff * p1) +
               (1 - p0) / (n0.eff * p0))

  list(
    est = logrr,
    se = se,
    low = logrr - 1.96 * se,
    up = logrr + 1.96 * se,
    p = 2 * (1 - pnorm(abs(logrr / se)))
  )
}

get.weighted.RR.sandwich <- function(data, w, gamma = NULL,
                                     target_mean = NULL,
                                     target_var = NULL,
                                     m_AD = NULL,
                                     target_cov = NULL,
                                     target_mean_vcov = NULL,
                                     account_for_gamma = TRUE,
                                     account_for_target_mean = FALSE) {
  x <- data$x
  y <- data$y
  w <- as.numeric(w)

  n <- length(y)

  idx1 <- x == 1
  idx0 <- x == 0
  I1 <- as.numeric(idx1)
  I0 <- as.numeric(idx0)

  W1 <- sum(w[idx1])
  W0 <- sum(w[idx0])

  p1 <- sum(w[idx1] * y[idx1]) / W1
  p0 <- sum(w[idx0] * y[idx0]) / W0

  use_gamma <- isTRUE(account_for_gamma)
  use_mu <- isTRUE(account_for_target_mean)

  if (use_mu && !use_gamma) {
    stop("To account for target_mean uncertainty, account_for_gamma should be TRUE.")
  }

  if (use_gamma) {
    if (is.null(target_mean)) {
      stop("target_mean must be supplied when account_for_gamma = TRUE.")
    }

    D <- sweep(covariate_matrix(data), 2, target_mean, "-")
    q <- ncol(D)

    Dbar.w <- colSums(w * D) / sum(w)
    Hgamma <- sweep(D, 2, Dbar.w, "-")

    psi.gamma <- w * D
    psi.p1 <- I1 * w * (y - p1)
    psi.p0 <- I0 * w * (y - p0)
    psi.ipd <- cbind(psi.gamma, psi.p1, psi.p0)

    A.gg <- crossprod(D, w * Hgamma) / n
    A.gp <- matrix(0, q, 2)
    A.pg <- rbind(
      colMeans(I1 * w * (y - p1) * Hgamma),
      colMeans(I0 * w * (y - p0) * Hgamma)
    )
    A.pp <- matrix(
      c(-mean(I1 * w), 0,
        0, -mean(I0 * w)),
      nrow = 2,
      byrow = TRUE
    )

    if (!use_mu) {
      A <- rbind(
        cbind(A.gg, A.gp),
        cbind(A.pg, A.pp)
      )

      B <- crossprod(psi.ipd) / n
      V <- MASS::ginv(A) %*% B %*% t(MASS::ginv(A)) / n

      grad <- c(rep(0, q), 1 / p1, -1 / p0)
      var.logrr <- as.numeric(t(grad) %*% V %*% grad)

      logrr <- log(p1 / p0)
      se <- sqrt(var.logrr)

      return(list(
        est = logrr,
        se = se,
        low = logrr - 1.96 * se,
        up = logrr + 1.96 * se,
        p = 2 * (1 - pnorm(abs(logrr / se)))
      ))
    }

    V.mu <- make_target_mean_vcov(
      target_var = target_var,
      m_AD = m_AD,
      target_cov = target_cov,
      target_mean_vcov = target_mean_vcov
    )

    if (!all(dim(V.mu) == c(q, q))) {
      stop("Dimension mismatch: Var(target_mean_hat) must be q by q.")
    }

    A.mu.mu <- -diag(q)
    A.mu.gamma <- matrix(0, q, q)
    A.mu.p <- matrix(0, q, 2)

    A.gamma.mu <- -mean(w) * diag(q)
    A.gamma.gamma <- A.gg
    A.gamma.p <- A.gp

    A.p.mu <- matrix(0, 2, q)
    A.p.gamma <- A.pg
    A.p.p <- A.pp

    A <- rbind(
      cbind(A.mu.mu,     A.mu.gamma,     A.mu.p),
      cbind(A.gamma.mu,  A.gamma.gamma,  A.gamma.p),
      cbind(A.p.mu,      A.p.gamma,      A.p.p)
    )

    B <- matrix(0, nrow = 2 * q + 2, ncol = 2 * q + 2)
    idx.mu <- seq_len(q)
    idx.ipd <- (q + 1):(2 * q + 2)

    B[idx.mu, idx.mu] <- n * V.mu
    B[idx.ipd, idx.ipd] <- crossprod(psi.ipd) / n

    V <- MASS::ginv(A) %*% B %*% t(MASS::ginv(A)) / n

    grad <- c(rep(0, q), rep(0, q), 1 / p1, -1 / p0)
    var.logrr <- as.numeric(t(grad) %*% V %*% grad)

    logrr <- log(p1 / p0)
    se <- sqrt(var.logrr)

    return(list(
      est = logrr,
      se = se,
      low = logrr - 1.96 * se,
      up = logrr + 1.96 * se,
      p = 2 * (1 - pnorm(abs(logrr / se))),
      V = V
    ))
  }

  psi <- cbind(
    I1 * w * (y - p1),
    I0 * w * (y - p0)
  )

  A <- matrix(
    c(-mean(I1 * w), 0,
      0, -mean(I0 * w)),
    nrow = 2,
    byrow = TRUE
  )

  B <- crossprod(psi) / n
  V <- MASS::ginv(A) %*% B %*% t(MASS::ginv(A)) / n

  grad <- c(1 / p1, -1 / p0)
  var.logrr <- as.numeric(t(grad) %*% V %*% grad)

  logrr <- log(p1 / p0)
  se <- sqrt(var.logrr)
  list(
    est = logrr,
    se = se,
    low = logrr - 1.96 * se,
    up = logrr + 1.96 * se,
    p = 2 * (1 - pnorm(abs(logrr / se)))
  )
}

bayes_est_RR_weighted <- function(data, w = NULL,
                                  a1 = 0.5, b1 = 0.5,
                                  a0 = 0.5, b0 = 0.5,
                                  M = 1e5, conf = 0.95) {
  if (is.null(w)) w <- rep(1, nrow(data))

  x <- data$x
  y <- data$y
  S1 <- sum(w[x == 1] * y[x == 1])
  F1 <- sum(w[x == 1] * (1 - y[x == 1]))
  S0 <- sum(w[x == 0] * y[x == 0])
  F0 <- sum(w[x == 0] * (1 - y[x == 0]))

  p1 <- rbeta(M, a1 + S1, b1 + F1)
  p0 <- rbeta(M, a0 + S0, b0 + F0)
  logrr <- log(p1 / p0)
  ci <- quantile(logrr, c((1 - conf) / 2, 1 - (1 - conf) / 2))

  list(
    point.est = mean(logrr),
    se.est = sd(logrr),
    conf.lower = unname(ci[1]),
    conf.upper = unname(ci[2]),
    p.value = NA_real_
  )
}

get.brm.ad <- function(est, data, weight) {
  if (is.null(est)) return(make_na_est())

  out <- est
  x <- data$x
  y <- data$y
  P1.w <- sum(weight[x == 1] * y[x == 1]) / sum(weight[x == 1])
  P0.w <- sum(weight[x == 0] * y[x == 0]) / sum(weight[x == 0])

  if (P0.w == 0 || P0.w == 1 || P1.w == 0 || P1.w == 1) {
    est.bayes <- bayes_est_RR_weighted(data, weight)
    out$point.est[1] <- est.bayes$point.est
    out$se.est[1] <- est.bayes$se.est
    out$conf.lower[1] <- est.bayes$conf.lower
    out$conf.upper[1] <- est.bayes$conf.upper
    out$p.value[1] <- NA_real_
  }

  out
}

gcomp_ml_rr_itc <- function(data.IPD, data.AD.or,
                            formula = NULL,
                            n_star = 500,
                            n_boot = 300,
                            correction = 0.5) {
  z.cols <- covariate_names(data.IPD$data)
  dat_ac <- data.IPD$data[, c("y", "x", z.cols), drop = FALSE]
  N.C <- data.AD.or$count[1]
  N.B <- data.AD.or$count[2]
  y.C.sum <- data.AD.or$count[3]
  y.B.sum <- data.AD.or$count[4]

  if (is.null(formula)) {
    main_terms <- paste(c("x", z.cols), collapse = " + ")
    interaction_terms <- paste(paste0("x:", z.cols), collapse = " + ")
    formula <- as.formula(paste("y ~", main_terms, "+", interaction_terms))
  }

  ad_z <- data.AD.or$data[, z.cols, drop = FALSE]
  x_star <- ad_z[sample.int(nrow(ad_z), n_star, replace = TRUE), , drop = FALSE]

  one_boot <- function(indices) {
    dat_b <- dat_ac[indices, , drop = FALSE]
    fit <- glm(formula, family = binomial(link = "logit"), data = dat_b)

    new_A <- x_star
    new_C <- x_star
    new_A$x <- 1
    new_C$x <- 0

    log(mean(predict(fit, newdata = new_A, type = "response")) /
          mean(predict(fit, newdata = new_C, type = "response")))
  }

  n <- nrow(dat_ac)
  boot_est <- replicate(n_boot, one_boot(sample.int(n, n, replace = TRUE)))

  pB <- y.B.sum / N.B
  pC <- y.C.sum / N.C
  if (any(c(y.B.sum, N.B - y.B.sum, y.C.sum, N.C - y.C.sum) == 0)) {
    y.B.sum <- y.B.sum + correction
    N.B <- N.B + 2 * correction
    y.C.sum <- y.C.sum + correction
    N.C <- N.C + 2 * correction
    pB <- y.B.sum / N.B
    pC <- y.C.sum / N.C
  }

  var.BC <- (1 / y.B.sum - 1 / N.B) + (1 / y.C.sum - 1 / N.C)
  est.AC <- mean(boot_est)
  var.AC <- var(boot_est)
  se.AC <- sqrt(var.AC)
  est.BC <- log(pB / pC)
  est.AB <- est.AC - est.BC
  se.AB <- sqrt(var.AC + var.BC)

  list(
    logRR.AC = est.AC,
    se.AC = se.AC,
    lower.AC = est.AC - 1.96 * se.AC,
    upper.AC = est.AC + 1.96 * se.AC,
    pvalAC = 2 * (1 - pnorm(abs(est.AC / se.AC))),
    logRR.BC = est.BC,
    logRR.AB = est.AB,
    se.AB = se.AB,
    lower.AB = est.AB - 1.96 * se.AB,
    upper.AB = est.AB + 1.96 * se.AB
  )
}

get.gcomp.ML <- function(data.IPD, data.AD.or,
                         n_star = 500,
                         n_boot = 300) {
  fit <- gcomp_ml_rr_itc(
    data.IPD = data.IPD,
    data.AD.or = data.AD.or,
    n_star = n_star,
    n_boot = n_boot
  )

  list(est = fit$logRR.AC, se = fit$se.AC,
       lower = fit$lower.AC, upper = fit$upper.AC, pval = fit$pvalAC)
}

exact_safe <- function(param, y, x, va, vb, weight,
                            point.est, se.est, pa, pb, label,
                            max.step = 40,
                            thres = 1e-3,
                            thres.dicho = 1e-3) {
  if (any(!is.finite(point.est[seq_len(pa + pb)]))) {
    warning(label, " skipped because point estimates are non-finite.",
            call. = FALSE)
    return(make_na_exact())
  }

  if (any(!is.finite(se.est[seq_len(pa)])) ||
      any(se.est[seq_len(pa)] <= 0)) {
    warning(label, " received non-finite alpha SE; using default exact-search width.",
            call. = FALSE)
    se.est[seq_len(pa)] <- 1
  }

  out <- safe_try(
    exact(param, y, x, va, vb, weight,
          max.step = max.step, thres = thres, thres.dicho = thres.dicho,
          pars = point.est, se = se.est, pa = pa, pb = pb),
    label
  )

  if (is.null(out)) make_na_exact() else out
}

get.compare.anchored <- function(data.ipd, data.ad, data.ad.or,
                                 param = "RR",
                                 run_exact = TRUE,
                                 include_gcomp = TRUE,
                                 n_star = 500,
                                 n_boot = 300) {
  y <- data.ipd$data$y
  x <- data.ipd$data$x
  va <- as.matrix(data.ipd$data$v.1, ncol = 1)
  vb <- cbind(data.ipd$data$v.1, covariate_matrix(data.ipd$data))
  pa <- ncol(va)
  pb <- ncol(vb)

  model.anweight <- get.weight(data.ipd, data.ad$anchored.mean.x)
  anchored.weight <- model.anweight$weight
  anchored.gamma <- model.anweight$ga

  brm.est <- safe_try(get.estimate("AD", data.ipd, rep(1, length(y)), 0, 0),
                      "brm")
  if (is.null(brm.est)) brm.est <- make_na_est()

  anchored.est <- safe_try(
    get.estimate("IPD", data.ipd, anchored.weight, anchored.gamma,
                 target_mean = data.ad$anchored.mean.x,
                 target_var = data.ad$anchored.var.x,
                 m_AD = data.ad$anchored.n,
                 target_cov = data.ad$anchored.cov.x),
    "brm_an"
  )
  if (is.null(anchored.est)) anchored.est <- make_na_est()

  LB.anchored <- safe_try(fit_logbinomial(data.ipd$data, anchored.weight,
                                          "LB.anchored"), "LB_an")
  LP.anchored <- safe_try(fit_lp(data.ipd$data, anchored.weight), "LP_an")
  RLP.anchored <- safe_try(fit_robust_lp(data.ipd$data, anchored.weight),
                           "RLP_an")
  CMH.anchored <- safe_try(get.weighted.RR(data.ipd$data, anchored.weight),
                           "CMH_an")
  CMH.sandwich.anchored <- safe_try(
    get.weighted.RR.sandwich(data.ipd$data, anchored.weight,
                             gamma = anchored.gamma,
                             target_mean = data.ad$anchored.mean.x,
                             target_var = data.ad$anchored.var.x,
                             m_AD = data.ad$anchored.n,
                             target_cov = data.ad$anchored.cov.x,
                             account_for_target_mean = TRUE),
    "CMH_sandwich_an"
  )

  bayes.anchored <- safe_try(get.brm.ad(anchored.est, data.ipd$data,
                                        anchored.weight), "brmad_an")
  if (is.null(bayes.anchored)) bayes.anchored <- make_na_est()

  est.gcomp <- if (isTRUE(include_gcomp)) {
    safe_try(get.gcomp.ML(data.ipd, data.ad.or, n_star, n_boot), "GC")
  } else {
    NULL
  }

  est.exact.anchored <- make_na_exact()
  est.exact.ad.anchored <- make_na_exact()
  if (isTRUE(run_exact)) {
    est.exact.anchored <- exact_safe(
      param, y, x, va, vb, anchored.weight,
      anchored.est$point.est, anchored.est$se.est, pa, pb, "exact_an"
    )
    est.exact.ad.anchored <- exact_safe(
      param, y, x, va, vb, anchored.weight,
      bayes.anchored$point.est, bayes.anchored$se.est, pa, pb, "ad_exact_an"
    )
  }

  result.comp <- rbind(
    point.est = c(
      stat_value(brm.est, "point.est"),
      stat_value(anchored.est, "point.est"),
      get_glm_stat(LB.anchored, "estimate"),
      get_glm_stat(LP.anchored, "estimate"),
      stat_value(RLP.anchored, "est"),
      stat_value(CMH.anchored, "est"),
      stat_value(CMH.sandwich.anchored, "est"),
      stat_value(bayes.anchored, "point.est"),
      stat_value(est.gcomp, "est"),
      stat_value(anchored.est, "point.est"),
      stat_value(bayes.anchored, "point.est")
    ),
    se.est = c(
      stat_value(brm.est, "se.est"),
      stat_value(anchored.est, "se.est"),
      get_glm_stat(LB.anchored, "se"),
      get_glm_stat(LP.anchored, "se"),
      stat_value(RLP.anchored, "se"),
      stat_value(CMH.anchored, "se"),
      stat_value(CMH.sandwich.anchored, "se"),
      stat_value(bayes.anchored, "se.est"),
      stat_value(est.gcomp, "se"),
      stat_value(anchored.est, "se.est"),
      stat_value(bayes.anchored, "se.est")
    ),
    con.low = c(
      stat_value(brm.est, "conf.lower"),
      stat_value(anchored.est, "conf.lower"),
      get_glm_stat(LB.anchored, "low"),
      get_glm_stat(LP.anchored, "low"),
      stat_value(RLP.anchored, "low"),
      stat_value(CMH.anchored, "low"),
      stat_value(CMH.sandwich.anchored, "low"),
      stat_value(bayes.anchored, "conf.lower"),
      stat_value(est.gcomp, "lower"),
      stat_value(est.exact.anchored, "low"),
      stat_value(est.exact.ad.anchored, "low")
    ),
    con.up = c(
      stat_value(brm.est, "conf.upper"),
      stat_value(anchored.est, "conf.upper"),
      get_glm_stat(LB.anchored, "up"),
      get_glm_stat(LP.anchored, "up"),
      stat_value(RLP.anchored, "up"),
      stat_value(CMH.anchored, "up"),
      stat_value(CMH.sandwich.anchored, "up"),
      stat_value(bayes.anchored, "conf.upper"),
      stat_value(est.gcomp, "upper"),
      stat_value(est.exact.anchored, "up"),
      stat_value(est.exact.ad.anchored, "up")
    ),
    p.value = c(
      stat_value(brm.est, "p.value"),
      stat_value(anchored.est, "p.value"),
      get_glm_stat(LB.anchored, "p"),
      get_glm_stat(LP.anchored, "p"),
      stat_value(RLP.anchored, "p"),
      stat_value(CMH.anchored, "p"),
      stat_value(CMH.sandwich.anchored, "p"),
      stat_value(bayes.anchored, "p.value"),
      stat_value(est.gcomp, "pval"),
      stat_value(est.exact.anchored, "p"),
      stat_value(est.exact.ad.anchored, "p")
    )
  )

  colnames(result.comp) <- c(
    "brm", "brm_an", "LB_an", "LP_an", "RLP_an",
    "CMH_an", "CMH_sandwich_an", "brmad_an", "GC", "exact_an", "ad_exact_an"
  )

  result.comp
}

run.anchored <- function(param, n, event, hypothesis,
                         run_exact = TRUE,
                         include_gcomp = TRUE,
                         n_star = 500,
                         n_boot = 300,
                         ess_ratio = NULL,
                         ess_side = "lower",
                         target_mean = NULL,
                         n_cov = 3,
                         ipd_prob = 0.5,
                         ad_prob = NULL,
                         direction = NULL) {
  data.IPD <- data.generation(
    param, "IPD", n, event, hypothesis,
    n_cov = n_cov,
    ipd_prob = ipd_prob
  )
  data.AD.or <- data.generation(
    param, "AD", n, event, hypothesis,
    ess_ratio = ess_ratio,
    ess_side = ess_side,
    target_mean = target_mean,
    n_cov = n_cov,
    ipd_prob = ipd_prob,
    ad_prob = ad_prob,
    direction = direction
  )
  data.AD <- get.aggre(data.AD.or)

  get.compare.anchored(
    data.IPD, data.AD, data.AD.or,
    param = param,
    run_exact = run_exact,
    include_gcomp = include_gcomp,
    n_star = n_star,
    n_boot = n_boot
  )
}
