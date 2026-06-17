## Clean unanchored helpers for ITC RR/RD simulation.
##
## Keep only the functions used by the final comparison:
##   weighted_un, LB_un, LP_un, RLP_un, GC_A_un, and optional brm_un.
##
## Required external files/functions before using brm_un/run.unanchored:
##   getProbScalarRR.R, MyFunc.R, 1.1_MLE_Point.R, 1.2_MLE_Var.R,
##   1.2_MLE_Var_joint.R, 1_CallMLE.R, and generate_RR.R.

## ------------------------------------------------------------------
## 0. Small utilities
## ------------------------------------------------------------------

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

clip_prob <- function(p, eps = 1e-8) {
  pmin(pmax(p, eps), 1 - eps)
}

if (!exists("log_sum_exp")) {
  log_sum_exp <- function(x) {
    x <- as.vector(x)
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
}

if (!exists("normalize_exp_weights")) {
  normalize_exp_weights <- function(eta) {
    eta <- as.vector(eta)
    w <- exp(eta - max(eta))
    w / mean(w)
  }
}

if (!exists("ginv")) {
  ginv <- function(X, ...) {
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::ginv(X, ...)
    } else {
      solve(X)
    }
  }
}

pinv <- function(A) ginv(A)

vcov_hc <- function(fit, type = "HC0") {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    sandwich::vcovHC(fit, type = type)
  } else {
    warning("Package 'sandwich' is unavailable; using model-based vcov().",
            call. = FALSE)
    stats::vcov(fit)
  }
}

covariate_names <- function(data) {
  grep("^z[0-9]+$", names(data), value = TRUE)
}

covariate_matrix <- function(data) {
  z.cols <- covariate_names(data)
  as.matrix(data[, z.cols, drop = FALSE])
}

ipd_arm_data <- function(data.IPD, arm = 1) {
  data.IPD$data[data.IPD$data$x == arm, , drop = FALSE]
}

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

## ------------------------------------------------------------------
## 1. Standard output format
## ------------------------------------------------------------------

UN_EST_ROWS <- c("point.est", "se.est", "con.low", "con.up", "p.value")
UN_EST_FIELDS <- c(
  point.est = "point.est",
  se.est    = "se.est",
  con.low   = "conf.lower",
  con.up    = "conf.upper",
  p.value   = "p.value"
)

make_na_est <- function() {
  list(
    point.est  = NA_real_,
    se.est     = NA_real_,
    conf.lower = NA_real_,
    conf.upper = NA_real_,
    p.value    = NA_real_
  )
}

make_est <- function(est, se) {
  est <- first_value(est)
  se <- first_value(se)

  if (!is.finite(est)) return(make_na_est())

  if (!is.finite(se) || se < 0) {
    return(list(
      point.est  = est,
      se.est     = NA_real_,
      conf.lower = NA_real_,
      conf.upper = NA_real_,
      p.value    = NA_real_
    ))
  }

  list(
    point.est  = est,
    se.est     = se,
    conf.lower = est - 1.96 * se,
    conf.upper = est + 1.96 * se,
    p.value    = if (se == 0) NA_real_ else 2 * (1 - pnorm(abs(est / se)))
  )
}

make_est_from_var <- function(est, var) {
  est <- first_value(est)
  var <- first_value(var)

  if (!is.finite(est)) return(make_na_est())
  if (!is.finite(var) || var < 0) {
    out <- make_est(est, NA_real_)
    out$var.est <- NA_real_
    return(out)
  }

  out <- make_est(est, sqrt(var))
  out$var.est <- var
  out
}

stat_value <- function(x, field) {
  if (is.null(x) || is.null(x[[field]])) return(NA_real_)
  first_value(x[[field]])
}

safe_est <- function(expr, label) {
  out <- safe_try(expr, label)
  if (is.null(out)) make_na_est() else out
}

est_to_column <- function(est) {
  out <- vapply(UN_EST_FIELDS, function(field) stat_value(est, field), numeric(1))
  names(out) <- UN_EST_ROWS
  out
}

est_matrix <- function(estimates) {
  out <- vapply(estimates, est_to_column, numeric(length(UN_EST_ROWS)))
  rownames(out) <- UN_EST_ROWS
  out
}

## ------------------------------------------------------------------
## 2. AD target summary and B-arm risk
## ------------------------------------------------------------------

binomial_risk <- function(N, y, correction = 0.5,
                          boundary = c("strict", "outer"),
                          clip = FALSE) {
  boundary <- match.arg(boundary)
  N.original <- N
  y.original <- y

  if (!is.finite(N) || N <= 0) {
    return(list(
      p = NA_real_, var = NA_real_, N = N, y = y,
      N.original = N.original, y.original = y.original
    ))
  }

  use.correction <- switch(
    boundary,
    strict = isTRUE(y == 0 || y == N),
    outer  = isTRUE(y <= 0 || y >= N)
  )

  if (use.correction) {
    y <- y + correction
    N <- N + 2 * correction
  }

  p <- y / N
  if (isTRUE(clip)) p <- clip_prob(p)

  list(
    p = p,
    var = p * (1 - p) / N,
    N = N,
    y = y,
    N.original = N.original,
    y.original = y.original
  )
}

get_AD_target_summary <- function(data.AD.or) {
  if (!is.null(data.AD.or$target.summary)) return(data.AD.or$target.summary)

  mean.x <- data.AD.or$target.mean.x %||%
    data.AD.or$unanchored.mean.x %||%
    data.AD.or$mean.x

  if (is.null(mean.x)) {
    stop(
      "AD target covariate summary is required for unanchored methods. ",
      "Provide data.AD.or$target.summary or target.mean.x/target.var.x/",
      "target.cov.x/target.n fields.",
      call. = FALSE
    )
  }

  cov.x <- data.AD.or$target.cov.x %||%
    data.AD.or$unanchored.cov.x %||%
    data.AD.or$cov.x

  var.x <- data.AD.or$target.var.x %||%
    data.AD.or$unanchored.var.x %||%
    data.AD.or$var.x

  n.AD <- data.AD.or$target.n %||%
    data.AD.or$unanchored.n %||%
    data.AD.or$n

  if (is.null(n.AD) && !is.null(data.AD.or$count)) n.AD <- data.AD.or$count[2]
  if (is.null(cov.x) && !is.null(var.x)) cov.x <- diag(as.numeric(var.x), length(mean.x))
  if (is.null(var.x) && !is.null(cov.x)) var.x <- diag(as.matrix(cov.x))

  list(
    mean.x = as.numeric(mean.x),
    var.x  = if (is.null(var.x)) NULL else as.numeric(var.x),
    cov.x  = if (is.null(cov.x)) NULL else as.matrix(cov.x),
    n      = n.AD
  )
}

target_mean_vcov <- function(target.summary) {
  q <- length(target.summary$mean.x)
  m.AD <- target.summary$n

  if (!is.null(target.summary$mean.vcov)) {
    return(as.matrix(target.summary$mean.vcov))
  }

  if (is.null(m.AD) || !is.finite(m.AD) || m.AD <= 1) {
    return(matrix(0, nrow = q, ncol = q))
  }

  if (!is.null(target.summary$cov.x) && all(is.finite(target.summary$cov.x))) {
    return(as.matrix(target.summary$cov.x) / m.AD)
  }

  if (!is.null(target.summary$var.x) && all(is.finite(target.summary$var.x))) {
    return(diag(as.numeric(target.summary$var.x) / m.AD, q))
  }

  matrix(0, nrow = q, ncol = q)
}

make_unanchored_ad_summary <- function(data.AD.or, data.IPD = NULL) {
  target.summary <- get_AD_target_summary(data.AD.or)
  q <- as.numeric(target.summary$mean.x)

  z.cols <- names(target.summary$mean.x)
  if (!is.null(data.IPD)) {
    z.ipd <- covariate_names(data.IPD$data)
    if (length(z.ipd) == length(q)) z.cols <- z.ipd
  }
  if (is.null(z.cols) || length(z.cols) != length(q)) {
    z.cols <- paste0("z", seq_along(q))
  }
  names(q) <- z.cols

  Vq <- target_mean_vcov(target.summary)
  if (!all(dim(Vq) == c(length(q), length(q)))) {
    stop("Dimension mismatch: Var(target mean) must be q by q.", call. = FALSE)
  }
  dimnames(Vq) <- list(z.cols, z.cols)

  N.B <- data.AD.or$N.B
  if (is.null(N.B) && !is.null(data.AD.or$count)) N.B <- data.AD.or$count[2]
  if (is.null(N.B)) N.B <- target.summary$N.B
  if (is.null(N.B)) N.B <- target.summary$n

  y.B.sum <- data.AD.or$y.B.sum
  if (is.null(y.B.sum) && !is.null(data.AD.or$count)) y.B.sum <- data.AD.or$count[4]
  if (is.null(y.B.sum)) y.B.sum <- target.summary$y.B.sum

  Cqp <- data.AD.or$C_z_pB %||% target.summary$C_z_pB %||% rep(0, length(q))
  Cqp <- matrix(as.numeric(Cqp), ncol = 1)
  if (nrow(Cqp) != length(q)) {
    stop("C_z_pB must have length equal to target mean length.", call. = FALSE)
  }

  list(
    target_mean    = q,
    var_target_mean = Vq,
    C_z_pB         = Cqp,
    N.B            = as.numeric(N.B),
    y.B.sum        = as.numeric(y.B.sum),
    z.cols         = z.cols
  )
}

get_B_info <- function(ad.sum, correction = 0.5) {
  risk <- binomial_risk(
    N = ad.sum$N.B,
    y = ad.sum$y.B.sum,
    correction = correction,
    boundary = "outer",
    clip = TRUE
  )

  list(pB = risk$p, Vp = risk$var, N.B = risk$N, y.B.sum = risk$y)
}

## ------------------------------------------------------------------
## 3. Common delta-method combination: pA versus pB
## ------------------------------------------------------------------

combine_unanchored_theory <- function(muA,
                                      V_mu_ipd,
                                      Gq,
                                      ad.sum,
                                      scale = c("logRR", "RD"),
                                      correction = 0.5,
                                      include_Cqp = TRUE) {
  scale <- match.arg(scale)
  muA <- clip_prob(muA)

  B.info <- get_B_info(ad.sum, correction = correction)
  pB <- B.info$pB
  Vp <- B.info$Vp

  Vq <- as.matrix(ad.sum$var_target_mean)
  Gq <- matrix(as.numeric(Gq), nrow = 1)

  if (ncol(Gq) == 0 || nrow(Vq) == 0) {
    Vq.part <- 0
    Cqp.part <- 0
  } else {
    Vq.part <- as.numeric(Gq %*% Vq %*% t(Gq))
    Cqp <- matrix(as.numeric(ad.sum$C_z_pB), ncol = 1)
    if (!isTRUE(include_Cqp)) Cqp[] <- 0
    Cqp.part <- as.numeric(Gq %*% Cqp)
  }

  if (scale == "logRR") {
    est <- log(muA) - log(pB)
    var.est <- V_mu_ipd / muA^2 + Vq.part / muA^2 +
      Vp / pB^2 - 2 * Cqp.part / (muA * pB)
  } else {
    est <- muA - pB
    var.est <- V_mu_ipd + Vq.part + Vp - 2 * Cqp.part
  }

  out <- make_est_from_var(est, var.est)
  out$muA <- muA
  out$pB <- pB
  out$V_mu_ipd <- as.numeric(V_mu_ipd)
  out$Vq_part <- as.numeric(Vq.part)
  out$Vp_part <- as.numeric(if (scale == "logRR") Vp / pB^2 else Vp)
  out$Cqp_part <- as.numeric(Cqp.part)
  out
}

## ------------------------------------------------------------------
## 4. MAIC-weighted one-arm estimator for pA
## ------------------------------------------------------------------

fit_maic_weights <- function(Z, target_mean,
                             warning_label = "Weight optimization did not fully converge: ") {
  Z <- as.matrix(Z)
  q <- ncol(Z)

  if (q == 0) {
    return(list(weight = rep(1, nrow(Z)), gamma = numeric(0), D = Z))
  }

  target_mean <- as.numeric(target_mean)
  if (length(target_mean) == 1 && q > 1) target_mean <- rep(target_mean, q)
  if (length(target_mean) != q) {
    stop("target_mean length must match the number of covariates.", call. = FALSE)
  }

  D <- sweep(Z, 2, target_mean, "-")

  Q <- function(alpha) log_sum_exp(D %*% alpha)
  Q.grad <- function(alpha) {
    eta <- as.vector(D %*% alpha)
    w <- exp(eta - max(eta))
    as.vector(crossprod(D, w / sum(w)))
  }

  fit <- optim(
    par = rep(0, q),
    fn = Q,
    gr = Q.grad,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  if (fit$convergence != 0) {
    warning(warning_label, fit$message, call. = FALSE)
  }

  list(
    weight = normalize_exp_weights(D %*% fit$par),
    gamma = fit$par,
    D = D,
    fit = fit
  )
}

fit_maic_weight_A <- function(datA, q, z.cols) {
  Z <- covariate_matrix(datA)[, z.cols, drop = FALSE]
  fit <- fit_maic_weights(
    Z,
    q,
    warning_label = "MAIC weight optimization may not have converged: "
  )
  list(weight = fit$weight, gamma = fit$gamma, D = fit$D)
}

estimate_weighted_un_theory <- function(data.IPD,
                                        ad.sum,
                                        scale = c("logRR", "RD"),
                                        correction = 0.5,
                                        include_Cqp = TRUE) {
  scale <- match.arg(scale)

  datA <- ipd_arm_data(data.IPD, 1)
  y <- datA$y
  nA <- nrow(datA)
  q <- ad.sum$target_mean
  z.cols <- ad.sum$z.cols

  wt <- fit_maic_weight_A(datA, q, z.cols)
  w <- wt$weight
  D <- wt$D
  mu <- clip_prob(sum(w * y) / sum(w))

  if (ncol(D) == 0) {
    psi <- matrix(w * (y - mu), ncol = 1)
    A <- matrix(-mean(w), 1, 1)
    B <- crossprod(psi) / nA
    Vtheta <- pinv(A) %*% B %*% t(pinv(A)) / nA
    V_mu_ipd <- as.numeric(Vtheta[1, 1])
    Gq <- numeric(0)
  } else {
    psi <- cbind(w * D, w * (y - mu))

    A11 <- crossprod(D, w * D) / nA
    A21 <- matrix(colMeans(w * (y - mu) * D), nrow = 1)
    A22 <- matrix(-mean(w), nrow = 1, ncol = 1)
    A <- rbind(
      cbind(A11, matrix(0, ncol(D), 1)),
      cbind(A21, A22)
    )

    B <- crossprod(psi) / nA
    Vtheta <- pinv(A) %*% B %*% t(pinv(A)) / nA
    V_mu_ipd <- as.numeric(Vtheta[ncol(Vtheta), ncol(Vtheta)])

    S_DD <- crossprod(D, w * D)
    S_YD <- matrix(colSums(w * (y - mu) * D), nrow = 1)
    Gq <- as.numeric(S_YD %*% pinv(S_DD))
  }

  combine_unanchored_theory(
    muA = mu,
    V_mu_ipd = V_mu_ipd,
    Gq = Gq,
    ad.sum = ad.sum,
    scale = scale,
    correction = correction,
    include_Cqp = include_Cqp
  )
}

## ------------------------------------------------------------------
## 5. One-arm GLM estimators for pA
## ------------------------------------------------------------------

target_row_from_mean <- function(target_mean, x = NULL, z.cols = names(target_mean)) {
  target_mean <- as.numeric(target_mean)
  if (is.null(z.cols) || length(z.cols) != length(target_mean)) {
    z.cols <- paste0("z", seq_along(target_mean))
  }

  dat <- as.data.frame(as.list(target_mean), check.names = FALSE)
  names(dat) <- z.cols
  dat$v.1 <- 1
  dat <- dat[, c("v.1", z.cols), drop = FALSE]
  if (!is.null(x)) dat$x <- x
  dat
}

one_arm_formula <- function(dat) {
  z.cols <- covariate_names(dat)
  if (length(z.cols) == 0) return(y ~ 1)
  as.formula(paste("y ~", paste(z.cols, collapse = " + ")))
}

one_arm_log_start <- function(form, dat, max_prob = 0.8, correction = 0.5) {
  mm <- model.matrix(form, data = dat)
  p <- (sum(dat$y) + correction) / (nrow(dat) + 2 * correction)
  p <- pmin(pmax(p, 1e-8), max_prob)

  start <- rep(0, ncol(mm))
  names(start) <- colnames(mm)
  start["(Intercept)"] <- log(p)
  start[is.na(start)] <- 0
  start
}

glm_family_for_method <- function(method = c("LB", "LP", "RLP", "LOGIT")) {
  method <- match.arg(method)
  switch(
    method,
    LB    = binomial(link = "log"),
    LP    = poisson(link = "log"),
    RLP   = poisson(link = "log"),
    LOGIT = binomial(link = "logit")
  )
}

fit_one_arm_glm <- function(datA, method = c("LB", "LP", "RLP", "LOGIT")) {
  method <- match.arg(method)
  form <- one_arm_formula(datA)
  fam <- glm_family_for_method(method)

  if (method == "LB") {
    glm(form, family = fam, data = datA,
        start = one_arm_log_start(form, datA, max_prob = 0.8))
  } else {
    glm(form, family = fam, data = datA)
  }
}

model_matrix_new <- function(fit, newdata) {
  tt <- stats::delete.response(stats::terms(fit))
  mf <- stats::model.frame(tt, newdata, xlev = fit$xlevels,
                           na.action = stats::na.pass)
  stats::model.matrix(tt, mf, contrasts.arg = fit$contrasts)
}

predict_mu_glm <- function(fit, newdata) {
  clip_prob(stats::predict(fit, newdata = newdata, type = "response"))
}

central_diff <- function(fn, par, eps = 1e-6) {
  par <- as.numeric(par)
  vapply(seq_along(par), function(j) {
    step <- eps * max(1, abs(par[j]))
    plus <- minus <- par
    plus[j] <- plus[j] + step
    minus[j] <- minus[j] - step
    (first_value(fn(plus)) - first_value(fn(minus))) / (2 * step)
  }, numeric(1))
}

grad_beta_glm <- function(fit, newdata) {
  Xnew <- model_matrix_new(fit, newdata)
  mu <- predict_mu_glm(fit, newdata)

  dmu.deta <- switch(
    fit$family$link,
    log      = mu,
    logit    = mu * (1 - mu),
    identity = rep(1, length(mu)),
    mu * (1 - mu)
  )

  as.numeric(dmu.deta) * Xnew
}

grad_q_glm_numeric <- function(fit, q, z.cols, x = NULL, eps = 1e-6) {
  d <- length(z.cols)
  if (d == 0) return(numeric(0))

  q <- as.numeric(q)
  names(q) <- z.cols

  central_diff(function(q.new) {
    names(q.new) <- z.cols
    predict_mu_glm(
      fit,
      target_row_from_mean(q.new, x = x, z.cols = z.cols)
    )
  }, q, eps = eps)
}

estimate_one_arm_glm_un_theory <- function(data.IPD,
                                           ad.sum,
                                           method = c("LB", "LP", "RLP", "LOGIT"),
                                           scale = c("logRR", "RD"),
                                           correction = 0.5,
                                           include_Cqp = TRUE,
                                           robust_vcov = NULL,
                                           hc_type = "HC0") {
  method <- match.arg(method)
  scale <- match.arg(scale)

  datA <- ipd_arm_data(data.IPD, 1)
  q <- ad.sum$target_mean
  z.cols <- ad.sum$z.cols

  fit <- fit_one_arm_glm(datA, method)
  newdat <- target_row_from_mean(q, z.cols = z.cols)
  muA <- first_value(predict_mu_glm(fit, newdat))

  Gbeta <- matrix(grad_beta_glm(fit, newdat), nrow = 1)
  if (is.null(robust_vcov)) robust_vcov <- method %in% c("RLP", "LOGIT")
  Vbeta <- if (isTRUE(robust_vcov)) vcov_hc(fit, type = hc_type) else vcov(fit)

  V_mu_ipd <- as.numeric(Gbeta %*% Vbeta %*% t(Gbeta))
  Gq <- grad_q_glm_numeric(fit, q, z.cols)

  out <- combine_unanchored_theory(
    muA = muA,
    V_mu_ipd = V_mu_ipd,
    Gq = Gq,
    ad.sum = ad.sum,
    scale = scale,
    correction = correction,
    include_Cqp = include_Cqp
  )

  out$fit <- fit
  out$Gq <- Gq
  out$V_mu_ipd <- V_mu_ipd
  out
}

## ------------------------------------------------------------------
## 6. Optional brm_un estimator based on weighted A + pseudo B data
## ------------------------------------------------------------------

get.estimate <- function(group, data, weight, gamma.ipd = NULL,
                         target_mean = 0,
                         target_var = NULL,
                         m_AD = NULL,
                         target_cov = NULL,
                         target_mean_vcov = NULL) {
  y <- data$data$y
  x <- data$data$x
  va <- matrix(data$data$v.1, ncol = 1)
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

as_count <- function(x, label) {
  x <- first_value(x)
  if (!is.finite(x) || x < 0) {
    stop(label, " must be a non-negative finite count.", call. = FALSE)
  }

  x.round <- round(x)
  if (abs(x - x.round) > 1e-8) {
    stop(label, " must be an integer count for pseudo-individual expansion.",
         call. = FALSE)
  }

  as.integer(x.round)
}

make_pseudo_B_data <- function(ad.sum) {
  N.B <- as_count(ad.sum$N.B, "N.B")
  y.B.sum <- as_count(ad.sum$y.B.sum, "y.B.sum")

  if (N.B <= 0) stop("N.B must be positive.", call. = FALSE)
  if (y.B.sum > N.B) stop("y.B.sum cannot exceed N.B.", call. = FALSE)

  data.frame(
    y = c(rep(1, y.B.sum), rep(0, N.B - y.B.sum)),
    x = rep(0, N.B),
    v.1 = rep(1, N.B),
    check.names = FALSE
  )
}

make_weighted_A_pseudo_B_data <- function(data.IPD, ad.sum) {
  datA <- ipd_arm_data(data.IPD, 1)
  if (nrow(datA) == 0) stop("No A-arm rows found in IPD data.", call. = FALSE)

  wt <- fit_maic_weight_A(datA, ad.sum$target_mean, ad.sum$z.cols)

  A.no.cov <- data.frame(
    y = datA$y,
    x = rep(1, nrow(datA)),
    v.1 = rep(1, nrow(datA)),
    check.names = FALSE
  )

  B.pseudo <- make_pseudo_B_data(ad.sum)

  list(
    data = list(data = rbind(A.no.cov, B.pseudo)),
    weight = c(wt$weight, rep(1, nrow(B.pseudo))),
    gamma = wt$gamma
  )
}

get.brm_un <- function(data.IPD, data.AD.or) {
  ad.sum <- make_unanchored_ad_summary(data.AD.or, data.IPD)
  dat <- make_weighted_A_pseudo_B_data(data.IPD, ad.sum)

  get.estimate(
    group = "AD",
    data = dat$data,
    weight = dat$weight,
    gamma.ipd = numeric(0),
    target_mean = numeric(0)
  )
}

brm_un_as_scale <- function(fit, scale = c("logRR", "RD")) {
  scale <- match.arg(scale)
  if (scale == "logRR") return(fit)

  theta <- stat_value(fit, "point.est")
  beta <- first_value(fit$point.est[2])
  V <- fit$cov

  if (!is.finite(theta) || !is.finite(beta) ||
      is.null(V) || !all(dim(V) == c(2, 2))) {
    return(make_na_est())
  }

  rd_fun <- function(par) {
    p <- getProbRR(par[1], par[2])
    first_value(p[, 2] - p[, 1])
  }

  par <- c(theta, beta)
  grad <- central_diff(rd_fun, par)
  make_est_from_var(rd_fun(par), as.numeric(t(grad) %*% V %*% grad))
}

## ------------------------------------------------------------------
## 7. Final comparison functions
## ------------------------------------------------------------------

get.compare.unanchored <- function(data.IPD,
                                   data.AD.or,
                                   scale = c("logRR", "RD"),
                                   include_brm_un = TRUE,
                                   correction = 0.5,
                                   include_Cqp = TRUE,
                                   hc_type = "HC0",
                                   ...) {
  scale <- match.arg(scale)
  ad.sum <- make_unanchored_ad_summary(data.AD.or, data.IPD)

  estimates <- list(
    weighted_un = safe_est(
      estimate_weighted_un_theory(
        data.IPD, ad.sum,
        scale = scale,
        correction = correction,
        include_Cqp = include_Cqp
      ),
      "weighted_un"
    ),
    LB_un = safe_est(
      estimate_one_arm_glm_un_theory(
        data.IPD, ad.sum,
        method = "LB",
        scale = scale,
        robust_vcov = FALSE,
        correction = correction,
        include_Cqp = include_Cqp,
        hc_type = hc_type
      ),
      "LB_un"
    ),
    LP_un = safe_est(
      estimate_one_arm_glm_un_theory(
        data.IPD, ad.sum,
        method = "LP",
        scale = scale,
        robust_vcov = FALSE,
        correction = correction,
        include_Cqp = include_Cqp,
        hc_type = hc_type
      ),
      "LP_un"
    ),
    RLP_un = safe_est(
      estimate_one_arm_glm_un_theory(
        data.IPD, ad.sum,
        method = "RLP",
        scale = scale,
        robust_vcov = TRUE,
        correction = correction,
        include_Cqp = include_Cqp,
        hc_type = hc_type
      ),
      "RLP_un"
    ),
    GC_A_un = safe_est(
      estimate_one_arm_glm_un_theory(
        data.IPD, ad.sum,
        method = "LOGIT",
        scale = scale,
        robust_vcov = TRUE,
        correction = correction,
        include_Cqp = include_Cqp,
        hc_type = hc_type
      ),
      "GC_A_un"
    )
  )

  if (isTRUE(include_brm_un)) {
    estimates$brm_un <- safe_est(
      brm_un_as_scale(get.brm_un(data.IPD, data.AD.or), scale = scale),
      "brm_un"
    )
  }

  est_matrix(estimates)
}

run.unanchored <- function(param,
                           n,
                           event,
                           hypothesis,
                           scale = c("logRR", "RD"),
                           include_brm_un = TRUE,
                           correction = 0.5,
                           include_Cqp = TRUE,
                           hc_type = "HC0",
                           ess_ratio = NULL,
                           ess_side = "lower",
                           target_mean = NULL,
                           n_cov = 3,
                           ipd_prob = 0.5,
                           ad_prob = NULL,
                           direction = NULL,
                           ...) {
  scale <- match.arg(scale)

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

  get.compare.unanchored(
    data.IPD = data.IPD,
    data.AD.or = data.AD.or,
    scale = scale,
    include_brm_un = include_brm_un,
    correction = correction,
    include_Cqp = include_Cqp,
    hc_type = hc_type
  )
}
