## ============================================================
## Unanchored ITC with AD summary data only
## AD supplies: target covariate mean q = E_B[Z], Var(qhat),
##              N_B and y_B_sum. Optional Cov(qhat, pBhat).
##
## Estimands:
##   logRR_un = log(mu_A^*) - log(p_B)
##   RD_un    = mu_A^* - p_B
##
## Methods included:
##   weighted_un      : MAIC-style weighted A-arm risk, sandwich for (gamma, mu)
##   LB_STC_un        : A-arm log-binomial STC at q
##   LP_STC_un        : A-arm Poisson log-link STC at q
##   RLP_STC_un       : A-arm robust Poisson STC at q
##   GC_A_STC_un      : A-arm logistic STC at q
##   GC_AC_STC_un     : A/C logistic STC at q, using IPD control as model assistance
##   brm_STC_un       : optional, needs MLEst/get.estimate to return alpha, beta and vcov
##
## Dependencies expected from your current project:
##   data.IPD$data must contain y, x, and covariates v.1, v.2, ...
##   x == 1 denotes IPD treatment A.
##
## Optional packages:
##   sandwich, MASS
## ============================================================

## --------------------------
## Basic helpers
## --------------------------

if (!exists("safe_try")) {
  safe_try <- function(expr, label) {
    tryCatch(expr, error = function(e) {
      warning(label, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }
}

if (!exists("make_na_est")) {
  make_na_est <- function() {
    list(point.est = NA_real_, se.est = NA_real_,
         conf.lower = NA_real_, conf.upper = NA_real_, p.value = NA_real_)
  }
}

if (!exists("first_value")) {
  first_value <- function(x) if (length(x) == 0) NA_real_ else x[1]
}

if (!exists("stat_value")) {
  stat_value <- function(x, field) {
    if (is.null(x) || is.null(x[[field]])) return(NA_real_)
    first_value(x[[field]])
  }
}

clip_prob <- function(p, eps = 1e-8) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

make_est <- function(est, var, scale = c("logRR", "RD")) {
  scale <- match.arg(scale)
  se <- sqrt(max(as.numeric(var), 0))
  list(
    point.est  = as.numeric(est),
    se.est     = se,
    conf.lower = as.numeric(est) - 1.96 * se,
    conf.upper = as.numeric(est) + 1.96 * se,
    p.value    = 2 * (1 - pnorm(abs(as.numeric(est) / se))),
    scale      = scale,
    var.est    = as.numeric(var)
  )
}

pinv <- function(A) {
  if (requireNamespace("MASS", quietly = TRUE)) return(MASS::ginv(A))
  solve(A)
}

vcov_hc <- function(fit, type = "HC0") {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    return(sandwich::vcovHC(fit, type = type))
  }
  warning("Package 'sandwich' is unavailable; using model-based vcov().", call. = FALSE)
  stats::vcov(fit)
}

## Prefer your existing covariate_names()/covariate_matrix() if they exist.
covariate_names_un <- function(dat) {
  if (exists("covariate_names")) return(covariate_names(dat))
  setdiff(names(dat), c("y", "x", "arm", "stratum"))
}

covariate_matrix_un <- function(dat, z.cols = NULL) {
  if (is.null(z.cols)) z.cols <- covariate_names_un(dat)
  if (length(z.cols) == 0) return(matrix(nrow = nrow(dat), ncol = 0))
  as.matrix(dat[, z.cols, drop = FALSE])
}

as_named_target_mean <- function(q, z.cols) {
  q <- as.numeric(q)
  if (length(q) != length(z.cols)) {
    stop("Length of target_mean does not match number of covariates.")
  }
  names(q) <- z.cols
  q
}

as_vcov_q <- function(Vq, d, N.B = NULL, var_is_for_mean = TRUE) {
  if (d == 0) return(matrix(0, 0, 0))

  if (is.null(Vq)) {
    stop("AD summary must supply Var(qhat), unless you set it to a zero matrix intentionally.")
  }

  if (is.vector(Vq) && !is.matrix(Vq)) {
    if (length(Vq) != d) stop("Length of var_target_mean vector is not d.")
    Vq <- diag(as.numeric(Vq), d)
  } else {
    Vq <- as.matrix(Vq)
    if (!all(dim(Vq) == c(d, d))) stop("var_target_mean must be d x d.")
  }

  if (!isTRUE(var_is_for_mean)) {
    if (is.null(N.B)) stop("N_B is needed if Vq is individual-level covariance.")
    Vq <- Vq / N.B
  }

  Vq
}

as_cov_q_p <- function(Cqp, d, N.B = NULL, cov_is_for_mean = TRUE) {
  if (d == 0) return(matrix(0, 0, 1))
  if (is.null(Cqp)) return(matrix(0, d, 1))

  Cqp <- matrix(as.numeric(Cqp), ncol = 1)
  if (nrow(Cqp) != d) stop("C_z_pB must have length d.")

  if (!isTRUE(cov_is_for_mean)) {
    if (is.null(N.B)) stop("N_B is needed if C_z_pB is individual-level covariance.")
    Cqp <- Cqp / N.B
  }

  Cqp
}

## --------------------------
## AD summary object
## --------------------------

make_AD_summary <- function(target_mean,
                            var_target_mean,
                            N_B,
                            y_B_sum,
                            z.cols = names(target_mean),
                            C_z_pB = NULL,
                            var_is_for_mean = TRUE,
                            cov_is_for_mean = TRUE) {
  if (is.null(z.cols)) z.cols <- paste0("v.", seq_along(target_mean))
  q <- as_named_target_mean(target_mean, z.cols)
  Vq <- as_vcov_q(var_target_mean, length(q), N.B = N_B,
                  var_is_for_mean = var_is_for_mean)
  Cqp <- as_cov_q_p(C_z_pB, length(q), N.B = N_B,
                    cov_is_for_mean = cov_is_for_mean)

  list(
    target_mean = q,
    var_target_mean = Vq,     # Var(qhat), not individual-level Var(Z), unless divided above
    C_z_pB = Cqp,             # Cov(qhat, pBhat); zero if not supplied
    N.B = as.numeric(N_B),
    y.B.sum = as.numeric(y_B_sum),
    z.cols = z.cols
  )
}

## Convenience function for simulations where data.AD.or still has individual AD data.
## In real use, call make_AD_summary() directly from published summary data.
make_AD_summary_from_AD_data <- function(data.AD.or,
                                         target_arm = 1,
                                         var_is_for_mean = TRUE) {
  dat <- data.AD.or$data
  z.cols <- covariate_names_un(dat)
  datB <- dat[dat$x == target_arm, , drop = FALSE]
  if (nrow(datB) == 0) stop("No B-arm rows found in AD individual-level simulation data.")

  ZB <- covariate_matrix_un(datB, z.cols)
  q <- colMeans(ZB)
  Vq <- stats::cov(ZB) / nrow(datB)
  if (length(z.cols) == 1) Vq <- matrix(stats::var(ZB[, 1]) / nrow(datB), 1, 1)

  ## In your current data structure: count = c(N.C, N.B, y.C.sum, y.B.sum).
  N.B <- data.AD.or$count[2]
  y.B.sum <- data.AD.or$count[4]

  make_AD_summary(q, Vq, N.B, y.B.sum, z.cols = z.cols,
                  C_z_pB = NULL, var_is_for_mean = TRUE)
}

get_B_info <- function(ad.sum, correction = 0.5) {
  N.B <- ad.sum$N.B
  yB <- ad.sum$y.B.sum

  if (yB <= 0 || yB >= N.B) {
    yB.tilde <- yB + correction
    N.tilde <- N.B + 2 * correction
  } else {
    yB.tilde <- yB
    N.tilde <- N.B
  }

  pB <- clip_prob(yB.tilde / N.tilde)
  Vp <- pB * (1 - pB) / N.tilde

  list(pB = pB, Vp = Vp, N.B = N.tilde, y.B.sum = yB.tilde)
}

## --------------------------
## Final variance combiner
## --------------------------

combine_unanchored_summary <- function(muA,
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
    var.est <- V_mu_ipd / muA^2 + Vq.part / muA^2 + Vp / pB^2 -
      2 * Cqp.part / (muA * pB)
  } else {
    est <- muA - pB
    var.est <- V_mu_ipd + Vq.part + Vp - 2 * Cqp.part
  }

  out <- make_est(est, var.est, scale = scale)
  out$muA <- muA
  out$pB <- pB
  out$V_mu_ipd <- as.numeric(V_mu_ipd)
  out$Vq_part <- as.numeric(Vq.part)
  out$Vp_part <- as.numeric(if (scale == "logRR") Vp / pB^2 else Vp)
  out$Cqp_part <- as.numeric(Cqp.part)
  out
}

## --------------------------
## 1. MAIC / weighted_un with sandwich and Var(qhat)
## --------------------------

fit_maic_weight_A <- function(datA, q, z.cols) {
  Z <- covariate_matrix_un(datA, z.cols)
  D <- sweep(Z, 2, q, "-")

  if (ncol(D) == 0) {
    return(list(weight = rep(1, nrow(datA)), gamma = numeric(0), D = D, Z = Z))
  }

  Q <- function(gamma) sum(exp(as.vector(D %*% gamma)))
  opt <- stats::optim(par = rep(0, ncol(D)), fn = Q, method = "BFGS")
  if (opt$convergence != 0) warning("MAIC weight optimization may not have converged.", call. = FALSE)

  w <- as.vector(exp(D %*% opt$par))
  w <- w / mean(w)
  list(weight = w, gamma = opt$par, D = D, Z = Z)
}

estimate_weighted_un <- function(data.IPD,
                                 ad.sum,
                                 scale = c("logRR", "RD"),
                                 correction = 0.5,
                                 include_Cqp = TRUE) {
  scale <- match.arg(scale)
  z.cols <- ad.sum$z.cols
  q <- ad.sum$target_mean

  datA <- data.IPD$data[data.IPD$data$x == 1, , drop = FALSE]
  y <- datA$y
  nA <- nrow(datA)

  wt <- fit_maic_weight_A(datA, q, z.cols)
  w <- wt$weight
  D <- wt$D

  mu <- sum(w * y) / sum(w)
  mu <- clip_prob(mu)

  if (ncol(D) == 0) {
    psi <- matrix(w * (y - mu), ncol = 1)
    A <- matrix(-mean(w), 1, 1)
    B <- crossprod(psi) / nA
    Vtheta <- pinv(A) %*% B %*% t(pinv(A)) / nA
    V_mu_ipd <- as.numeric(Vtheta[1, 1])
    Gq <- numeric(0)
  } else {
    ## Joint estimating equation:
    ##   mean{w(gamma)(Z-q)} = 0
    ##   mean{w(gamma)(Y-mu)} = 0
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

    ## Derivative of mu_A^*(q) w.r.t. q, accounting for gamma(q).
    S_DD <- crossprod(D, w * D)
    S_YD <- matrix(colSums(w * (y - mu) * D), nrow = 1)
    Gq <- as.numeric(S_YD %*% pinv(S_DD))
  }

  out <- combine_unanchored_summary(
    muA = mu,
    V_mu_ipd = V_mu_ipd,
    Gq = Gq,
    ad.sum = ad.sum,
    scale = scale,
    correction = correction,
    include_Cqp = include_Cqp
  )

  out$method <- "weighted_un"
  out$gamma <- wt$gamma
  out
}

## --------------------------
## 2. Summary-only STC / GLM methods
## --------------------------

make_target_row <- function(q, z.cols, x_value = NULL) {
  if (length(z.cols) == 0) {
    newdat <- data.frame(.dummy = 1)
    newdat$.dummy <- NULL
  } else {
    newdat <- as.data.frame(as.list(q[z.cols]))
    names(newdat) <- z.cols
  }
  if (!is.null(x_value)) newdat$x <- x_value
  newdat
}

one_arm_formula <- function(z.cols) {
  if (length(z.cols) == 0) return(stats::as.formula("y ~ 1"))
  stats::as.formula(paste("y ~", paste(z.cols, collapse = " + ")))
}

ac_formula <- function(z.cols, include_interaction = TRUE) {
  if (length(z.cols) == 0) return(stats::as.formula("y ~ x"))

  main_terms <- paste(c("x", z.cols), collapse = " + ")
  if (isTRUE(include_interaction)) {
    int_terms <- paste(paste0("x:", z.cols), collapse = " + ")
    stats::as.formula(paste("y ~", main_terms, "+", int_terms))
  } else {
    stats::as.formula(paste("y ~", main_terms))
  }
}

model_matrix_new <- function(fit, newdata) {
  tt <- stats::delete.response(stats::terms(fit))
  mf <- stats::model.frame(tt, newdata, xlev = fit$xlevels, na.action = stats::na.pass)
  stats::model.matrix(tt, mf, contrasts.arg = fit$contrasts)
}

predict_mu_glm <- function(fit, newdata) {
  clip_prob(stats::predict(fit, newdata = newdata, type = "response"))
}

grad_beta_glm <- function(fit, newdata) {
  Xnew <- model_matrix_new(fit, newdata)
  mu <- predict_mu_glm(fit, newdata)
  link <- fit$family$link

  dmu_deta <- switch(
    link,
    log = mu,
    logit = mu * (1 - mu),
    identity = rep(1, length(mu)),
    ## fallback numerical derivative if another link appears
    mu * (1 - mu)
  )

  as.numeric(dmu_deta) * Xnew
}

grad_q_glm_numeric <- function(fit, q, z.cols, x_value = NULL,
                               eps = 1e-6) {
  d <- length(z.cols)
  if (d == 0) return(numeric(0))

  q <- as_named_target_mean(q, z.cols)
  g <- numeric(d)

  for (j in seq_len(d)) {
    step <- eps * max(1, abs(q[j]))
    q_plus <- q
    q_minus <- q
    q_plus[j] <- q_plus[j] + step
    q_minus[j] <- q_minus[j] - step

    nd_plus <- make_target_row(q_plus, z.cols, x_value = x_value)
    nd_minus <- make_target_row(q_minus, z.cols, x_value = x_value)

    mu_plus <- predict_mu_glm(fit, nd_plus)
    mu_minus <- predict_mu_glm(fit, nd_minus)
    g[j] <- (mu_plus - mu_minus) / (2 * step)
  }

  g
}

fit_stc_glm <- function(data.IPD,
                        ad.sum,
                        method = c("LB", "LP", "RLP", "LOGIT", "AC_LOGIT"),
                        scale = c("logRR", "RD"),
                        robust_vcov = NULL,
                        hc_type = "HC0",
                        correction = 0.5,
                        include_Cqp = TRUE,
                        include_interaction = TRUE) {
  method <- match.arg(method)
  scale <- match.arg(scale)

  z.cols <- ad.sum$z.cols
  q <- ad.sum$target_mean

  if (method == "AC_LOGIT") {
    dat <- data.IPD$data[, c("y", "x", z.cols), drop = FALSE]
    form <- ac_formula(z.cols, include_interaction = include_interaction)
    newdat <- make_target_row(q, z.cols, x_value = 1)
    fam <- stats::binomial(link = "logit")
    method_name <- "GC_AC_STC_un"
    if (is.null(robust_vcov)) robust_vcov <- TRUE
  } else {
    dat <- data.IPD$data[data.IPD$data$x == 1, c("y", z.cols), drop = FALSE]
    form <- one_arm_formula(z.cols)
    newdat <- make_target_row(q, z.cols, x_value = NULL)

    fam <- switch(
      method,
      LB    = stats::binomial(link = "log"),
      LP    = stats::poisson(link = "log"),
      RLP   = stats::poisson(link = "log"),
      LOGIT = stats::binomial(link = "logit")
    )

    method_name <- switch(
      method,
      LB = "LB_STC_un",
      LP = "LP_STC_un",
      RLP = "RLP_STC_un",
      LOGIT = "GC_A_STC_un"
    )

    if (is.null(robust_vcov)) robust_vcov <- method %in% c("RLP", "LOGIT")
  }

  fit <- if (method == "LB") {
    ## A conservative start helps log-binomial convergence.
    X0 <- stats::model.matrix(form, data = dat)
    stats::glm(form, family = fam, data = dat, start = rep(-0.01, ncol(X0)))
  } else {
    stats::glm(form, family = fam, data = dat)
  }

  muA <- predict_mu_glm(fit, newdat)
  Gbeta <- matrix(grad_beta_glm(fit, newdat), nrow = 1)

  Vbeta <- if (isTRUE(robust_vcov)) vcov_hc(fit, type = hc_type) else stats::vcov(fit)
  V_mu_ipd <- as.numeric(Gbeta %*% Vbeta %*% t(Gbeta))

  Gq <- grad_q_glm_numeric(
    fit = fit,
    q = q,
    z.cols = z.cols,
    x_value = if (method == "AC_LOGIT") 1 else NULL
  )

  out <- combine_unanchored_summary(
    muA = muA,
    V_mu_ipd = V_mu_ipd,
    Gq = Gq,
    ad.sum = ad.sum,
    scale = scale,
    correction = correction,
    include_Cqp = include_Cqp
  )

  out$method <- method_name
  out$fit <- fit
  out$Gq <- Gq
  out$V_mu_ipd <- V_mu_ipd
  out
}

## --------------------------
## 3. Optional brm_STC_un
## --------------------------
## This requires your MLEst()/get.estimate() to return:
##   alpha.ml or alpha or alpha.hat
##   beta.ml  or beta  or beta.hat
##   and a full vcov matrix for c(alpha, beta): vcov, var, var.mat, V, or cov.
## If your current MLEst only returns point.est and se.est, modify it first.

`%||%` <- function(a, b) if (!is.null(a)) a else b

extract_brm_par_vcov <- function(fit) {
  alpha <- fit$alpha.ml %||% fit$alpha %||% fit$alpha.hat
  beta  <- fit$beta.ml  %||% fit$beta  %||% fit$beta.hat
  Veta  <- fit$vcov %||% fit$var %||% fit$var.mat %||% fit$V %||% fit$cov

  if (is.null(alpha) || is.null(beta) || is.null(Veta)) {
    stop("brm_STC_un needs alpha, beta and full vcov for c(alpha, beta).")
  }

  list(alpha = as.numeric(alpha), beta = as.numeric(beta), Veta = as.matrix(Veta))
}

brm_design_q <- function(q, z.cols, effect_col = "v.1") {
  if (!(effect_col %in% z.cols)) effect_col <- z.cols[1]

  ## Mirrors your current get.estimate():
  ##   va <- as.matrix(data$v.1, ncol = 1)
  ##   vb <- cbind(data$v.1, covariate_matrix(data))
  va <- matrix(q[effect_col], nrow = 1)
  vb <- matrix(c(q[effect_col], q[z.cols]), nrow = 1)

  J_va <- matrix(0, nrow = ncol(va), ncol = length(z.cols))
  J_va[1, match(effect_col, z.cols)] <- 1

  J_vb <- matrix(0, nrow = ncol(vb), ncol = length(z.cols))
  J_vb[1, match(effect_col, z.cols)] <- 1
  J_vb[seq_along(z.cols) + 1, seq_along(z.cols)] <- diag(length(z.cols))

  list(va = va, vb = vb, J_va = J_va, J_vb = J_vb)
}

brm_rr_p1_grad <- function(theta, phi) {
  r <- exp(theta)
  o <- exp(phi)

  if (abs(o - 1) < 1e-8) {
    p0 <- 1 / (1 + r)
  } else {
    ## Solve: o(1-p0)(1-r*p0) - r*p0^2 = 0
    a <- r * (o - 1)
    b <- -o * (1 + r)
    c <- o
    disc <- max(b^2 - 4 * a * c, 0)
    root1 <- (-b - sqrt(disc)) / (2 * a)
    root2 <- (-b + sqrt(disc)) / (2 * a)
    upper <- min(1, 1 / r)
    p0 <- if (root1 > 0 && root1 < upper) root1 else root2
  }

  p0 <- clip_prob(p0)
  p1 <- clip_prob(r * p0)

  H <- 2 * r * p0 + o * (1 + r - 2 * r * p0)
  dp0_dtheta <- -r * (p0^2 + o * p0 * (1 - p0)) / H
  dp0_dphi <- o * (1 - p0) * (1 - r * p0) / H

  dp1_dtheta <- p1 + r * dp0_dtheta
  dp1_dphi <- r * dp0_dphi

  list(p1 = p1, dp1_dtheta = as.numeric(dp1_dtheta), dp1_dphi = as.numeric(dp1_dphi))
}

estimate_brm_stc_un <- function(data.IPD,
                                ad.sum,
                                scale = c("logRR", "RD"),
                                correction = 0.5,
                                include_Cqp = TRUE,
                                effect_col = "v.1") {
  scale <- match.arg(scale)

  if (!exists("get.estimate")) {
    stop("get.estimate() is not found. Source your BRM code first.")
  }

  y <- data.IPD$data$y
  fit <- get.estimate("AD", data.IPD, rep(1, length(y)), gamma.ipd = 0, target_mean = 0)
  pb <- extract_brm_par_vcov(fit)

  z.cols <- ad.sum$z.cols
  q <- ad.sum$target_mean
  des <- brm_design_q(q, z.cols, effect_col = effect_col)

  theta <- as.numeric(des$va %*% pb$alpha)
  phi <- as.numeric(des$vb %*% pb$beta)
  inv <- brm_rr_p1_grad(theta, phi)
  muA <- inv$p1

  G_eta <- c(inv$dp1_dtheta * as.numeric(des$va),
             inv$dp1_dphi   * as.numeric(des$vb))
  G_eta <- matrix(G_eta, nrow = 1)

  if (!all(dim(pb$Veta) == c(ncol(G_eta), ncol(G_eta)))) {
    stop("Dimension mismatch: brm vcov must match length c(alpha, beta).")
  }

  V_mu_ipd <- as.numeric(G_eta %*% pb$Veta %*% t(G_eta))

  dtheta_dq <- as.numeric(t(pb$alpha) %*% des$J_va)
  dphi_dq <- as.numeric(t(pb$beta) %*% des$J_vb)
  Gq <- inv$dp1_dtheta * dtheta_dq + inv$dp1_dphi * dphi_dq

  out <- combine_unanchored_summary(
    muA = muA,
    V_mu_ipd = V_mu_ipd,
    Gq = Gq,
    ad.sum = ad.sum,
    scale = scale,
    correction = correction,
    include_Cqp = include_Cqp
  )

  out$method <- "brm_STC_un"
  out$fit <- fit
  out$Gq <- Gq
  out$V_mu_ipd <- V_mu_ipd
  out
}

## --------------------------
## 4. Main comparison function
## --------------------------

get.compare.unanchored.summary <- function(data.IPD,
                                           ad.sum,
                                           scale = c("logRR", "RD"),
                                           correction = 0.5,
                                           include_Cqp = TRUE,
                                           include_brm = FALSE,
                                           hc_type = "HC0") {
  scale <- match.arg(scale)

  methods <- list()

  methods$weighted_un <- safe_try(
    estimate_weighted_un(data.IPD, ad.sum, scale = scale,
                         correction = correction, include_Cqp = include_Cqp),
    "weighted_un"
  )

  methods$LB_STC_un <- safe_try(
    fit_stc_glm(data.IPD, ad.sum, method = "LB", scale = scale,
                robust_vcov = FALSE, hc_type = hc_type,
                correction = correction, include_Cqp = include_Cqp),
    "LB_STC_un"
  )

  methods$LP_STC_un <- safe_try(
    fit_stc_glm(data.IPD, ad.sum, method = "LP", scale = scale,
                robust_vcov = FALSE, hc_type = hc_type,
                correction = correction, include_Cqp = include_Cqp),
    "LP_STC_un"
  )

  methods$RLP_STC_un <- safe_try(
    fit_stc_glm(data.IPD, ad.sum, method = "RLP", scale = scale,
                robust_vcov = TRUE, hc_type = hc_type,
                correction = correction, include_Cqp = include_Cqp),
    "RLP_STC_un"
  )

  methods$GC_A_STC_un <- safe_try(
    fit_stc_glm(data.IPD, ad.sum, method = "LOGIT", scale = scale,
                robust_vcov = TRUE, hc_type = hc_type,
                correction = correction, include_Cqp = include_Cqp),
    "GC_A_STC_un"
  )

  methods$GC_AC_STC_un <- safe_try(
    fit_stc_glm(data.IPD, ad.sum, method = "AC_LOGIT", scale = scale,
                robust_vcov = TRUE, hc_type = hc_type,
                correction = correction, include_Cqp = include_Cqp,
                include_interaction = TRUE),
    "GC_AC_STC_un"
  )

  if (isTRUE(include_brm)) {
    methods$brm_STC_un <- safe_try(
      estimate_brm_stc_un(data.IPD, ad.sum, scale = scale,
                          correction = correction, include_Cqp = include_Cqp),
      "brm_STC_un"
    )
  }

  methods <- lapply(methods, function(x) if (is.null(x)) make_na_est() else x)

  result.comp <- rbind(
    point.est = sapply(methods, stat_value, field = "point.est"),
    se.est    = sapply(methods, stat_value, field = "se.est"),
    con.low   = sapply(methods, stat_value, field = "conf.lower"),
    con.up    = sapply(methods, stat_value, field = "conf.upper"),
    p.value   = sapply(methods, stat_value, field = "p.value")
  )

  result.comp
}

## --------------------------
## 5. Example wrapper for simulations
## --------------------------
## Use this only if your simulated AD object still has individual data.
## In a real AD-summary-only analysis, create ad.sum by make_AD_summary().

run.unanchored.summary <- function(param, n, event, hypothesis,
                                   scale = c("logRR", "RD"),
                                   include_brm = FALSE,
                                   ess_ratio = NULL,
                                   ess_side = "lower",
                                   target_mean = NULL,
                                   n_cov = 3,
                                   ipd_prob = 0.5,
                                   ad_prob = NULL,
                                   direction = NULL) {
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

  ad.sum <- make_AD_summary_from_AD_data(data.AD.or, target_arm = 1)

  get.compare.unanchored.summary(
    data.IPD = data.IPD,
    ad.sum = ad.sum,
    scale = scale,
    include_brm = include_brm
  )
}
