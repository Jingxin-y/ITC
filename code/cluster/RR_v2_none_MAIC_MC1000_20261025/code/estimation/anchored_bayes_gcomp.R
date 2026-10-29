## 边界事件的贝叶斯修正，以及三种 BRM 的 BC 区间。
## 保留原文件名；已移除 G-computation。
## 方法组合参考 brmplus/compare/run_brmfirth.R。


get.estimate <- function(group, data, weight, gamma.ipd = NULL,
                         target_mean = 0,
                         target_var = NULL,
                         m_AD = NULL,
                         target_cov = NULL,
                         target_mean_vcov = NULL, method = "brm",
                         effect_modifiers = NULL, param = "RR",
                         weight_covariates = NULL, model_covariates = NULL,
                         non_effect_modifiers = NULL,
                         op_model = c("marginal", "linear"), score_spec = NULL,
                         fit_control = list()) {
  fit_control <- anchored_fit_control(fit_control)
  prepared <- prepare_brm_score_model(data, score_spec, effect_modifiers,
                                      non_effect_modifiers, model_covariates)
  data <- prepared$data
  score_spec <- prepared$score_spec
  effect_modifiers <- prepared$effect_modifiers
  non_effect_modifiers <- prepared$non_effect_modifiers
  model_covariates <- prepared$model_covariates
  y <- data$data$y
  x <- data$data$x
  ## 默认：RR/RD 保留线性 EM，OP 对同一 EM 的联合格饱和，非 EM 被 marginalize。
  ## OR 的旧 logistic 接口不使用 OP 参数化，继续沿用原线性设计。
  op_model <- if (param == "OR") "linear" else match.arg(op_model)
  design <- brm_op_design(data$data, effect_modifiers, non_effect_modifiers,
                          model_covariates, op_model)
  va <- design$va
  vb <- design$vb
  if (length(weight) != length(y) || any(!is.finite(weight)) ||
      any(weight < 0) || !any(weight > 0))
    stop("weight must contain finite nonnegative values with positive total weight.")
  if (anyNA(x) || any(!x %in% c(0, 1))) stop("x must be coded as 0/1.")
  if (op_model == "marginal") {
    k <- ncol(vb)
    n0 <- tabulate(design$cell_id[weight > 0 & x == 0], nbins = k)
    n1 <- tabulate(design$cell_id[weight > 0 & x == 1], nbins = k)
    design$spec$cells$n_arm0 <- n0
    design$spec$cells$n_arm1 <- n1
    informative <- weight > 0 & design$cell_id %in% which(n0 > 0 & n1 > 0)
    if (!any(informative) || qr(va[informative, , drop = FALSE])$rank < ncol(va))
      stop("Saturated OP needs enough EM cells with both treatment arms to identify the effect model; use discrete EM or a prespecified grouped score.")
  }
  if (qr(va[weight > 0, , drop = FALSE])$rank < ncol(va) ||
      qr(vb[weight > 0, , drop = FALSE])$rank < ncol(vb))
    stop("BRM design is not identifiable on positive-weight observations; check EM variation and OP cells.")
  
  pa <- ncol(va)
  pb <- ncol(vb)
  weight_covariates <- weight_covariate_names(data$data, weight_covariates, gamma.ipd)
  vg <- covariate_matrix(data$data, weight_covariates)
  q <- ncol(vg)
  all_covariates <- covariate_names(data$data)
  target_mean <- select_weight_summary(target_mean, weight_covariates, all_covariates)
  target_var <- select_weight_summary(target_var, weight_covariates, all_covariates)
  target_cov <- select_weight_summary(target_cov, weight_covariates, all_covariates)
  target_mean_vcov <- select_weight_summary(target_mean_vcov, weight_covariates, all_covariates)
  if (is.null(gamma.ipd)) gamma.ipd <- rep(0, q)
  if (!is.null(names(gamma.ipd))) gamma.ipd <- gamma.ipd[weight_covariates]
  if (length(gamma.ipd) == 1 && q > 1) gamma.ipd <- rep(gamma.ipd, q)
  
  fit <- MLEst(
    group, param, y, x, va, vb, weight, gamma.ipd,
    max.step = if (method == "firth") fit_control$fc_max_step else fit_control$brm_max_step,
    thres = if (method == "firth") fit_control$fc_thres else fit_control$brm_thres,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb),
    pa = pa,
    pb = pb,
    target_mean = target_mean,
    target_var = target_var,
    m_AD = m_AD,
    target_cov = target_cov,
    target_mean_vcov = target_mean_vcov, method = method, vg = vg,
    optim.maxit = fit_control$brm_optim_maxit,
    optim.reltol = fit_control$brm_optim_reltol
  )
  fit$op_spec <- design$spec
  fit$score_spec <- score_spec
  fit$fit_control <- fit_control
  fit
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
    if (!is.null(est$va) && ncol(est$va) > 1L) {
      stop("The arm-level Bayes fallback does not estimate effect-modifier coefficients.")
    }
    est.bayes <- bayes_est_RR_weighted(data, weight)
    out$point.est[1] <- est.bayes$point.est
    out$se.est[1] <- est.bayes$se.est
    out$conf.lower[1] <- est.bayes$conf.lower
    out$conf.upper[1] <- est.bayes$conf.upper
    out$p.value[1] <- NA_real_
    if (!is.null(out$coefficients)) {
      out$coefficients[1, ] <- c(est.bayes$point.est, est.bayes$se.est,
                                 est.bayes$conf.lower, est.bayes$conf.upper, NA_real_)
    }
    if (!is.null(out$cov)) {
      out$cov[1, ] <- out$cov[, 1] <- NA_real_
      out$cov[1, 1] <- est.bayes$se.est^2
    }
    if (!is.null(out$va)) out$param.est <- exp(out$va %*% out$point.est[1])
    out$bayes_adjusted <- TRUE
    out$covariance_note <- "Alpha-beta covariance is unavailable after the arm-level Bayes correction."
  }

  out
}

get.brm.bc <- function(data.ipd, data.ad, weight, gamma,
                       brm.est, bayes.est, run_exact = TRUE,
                       seed = NULL, exact_seed_offset = 1000000L,
                       bc_control = list(), weight_covariates = NULL,
                       model_covariates = if (!is.null(brm.est$op_spec))
                         brm.est$op_spec$model_covariates else colnames(brm.est$vb)[-1],
                       effect_modifiers = colnames(brm.est$va)[-1],
                       non_effect_modifiers = brm.est$op_spec$non_effect_modifiers,
                       op_model = if (!is.null(brm.est$op_spec)) brm.est$op_spec$mode else "linear",
                       score_spec = brm.est$score_spec, fit_control = list(),
                       bc_target = c("marginal", "coefficient"),
                       model_weighting = c("unit", "maic"),
                       standardization.weight = weight) {
  model_weighting <- match.arg(model_weighting)
  ## Capture MAIC target weights before replacing likelihood weights with ones.
  ## Standalone calls default to the same unit-fit / MAIC-target convention as
  ## the main pipeline; weighted fits require explicit model_weighting="maic".
  standardization.weight <- as.numeric(standardization.weight)
  if (model_weighting == "unit") {
    weight <- rep(1, nrow(data.ipd$data))
    gamma <- setNames(numeric(0), character(0))
    weight_covariates <- character(0)
  }
  model.group <- if (model_weighting == "unit") "AD" else "IPD"
  bc_target <- match.arg(bc_target)
  if (bc_target == "marginal") validate_marginal_bc_control(bc_control)
  ## 标量 Beta-Binomial 修正不用于多参数治疗效应模型。
  if (length(effect_modifiers)) bayes.est <- make_na_est()
  y <- data.ipd$data$y
  x <- data.ipd$data$x
  if (!is.null(score_spec)) {
    score_weights <- intersect(weight_covariate_names(data.ipd$data, weight_covariates, gamma),
                                c("S_E","S_P","h_E","h_P"))
    data.ad <- brm_score_summary(data.ad, score_spec, score_weights)
  }

  firth <- safe_try(
    get.estimate(model.group, data.ipd, weight, gamma,
                 target_mean = data.ad$anchored.mean.x,
                 target_var = data.ad$anchored.var.x,
                 m_AD = data.ad$anchored.n,
                 target_cov = data.ad$anchored.cov.x, method = "firth",
                 weight_covariates = weight_covariates,
                 model_covariates = model_covariates, effect_modifiers = effect_modifiers,
                 non_effect_modifiers = non_effect_modifiers, op_model = op_model,
                 score_spec = score_spec, fit_control = fit_control),
    "brm-FC"
  )
  if (is.null(firth)) firth <- make_na_est()

  fits <- list(brm.bc = brm.est, firth.bc = firth, bayes.bc = bayes.est)
  labels <- c("brm-BC", "brm-FC-BC", "brm_b-BC")
  if (bc_target == "marginal") {
    if (length(standardization.weight) != length(y) ||
        any(!is.finite(standardization.weight)) || any(standardization.weight < 0) ||
        !any(standardization.weight > 0)) {
      intervals <- lapply(fits, function(fit) make_na_marginal_exact("weighting_failed"))
      return(c(list(firth=firth), intervals))
    }
    intervals <- lapply(seq_along(fits), function(i) marginal_exact_from_fit(
      fits[[i]], y, x, weight=weight, standardization.weight=standardization.weight,
      run_exact=run_exact, seed=if(is.null(seed)) NULL else seed+exact_seed_offset,
      bc_control=bc_control, label=labels[i]))
    names(intervals) <- names(fits)
    return(c(list(firth=firth),intervals))
  }
  intervals <- lapply(fits, function(fit) make_na_exact(1L + length(effect_modifiers)))
  if (run_exact) {
    for (i in seq_along(fits)) {
      ## 与上游一致：每种方法的 BC 从同一个 seed + offset 开始。
      if (!is.null(seed)) set.seed(seed + exact_seed_offset)
      fit <- fits[[i]]
       if (is.null(fit$va) || is.null(fit$vb)) next # 拟合失败时保留 NA 区间。
       if (!isTRUE(fit$converged)) next
      ## BC 必须沿用该次拟合的完整设计矩阵。
      intervals[[i]] <- exact_safe(
        "RR", y, x, fit$va, fit$vb, weight, fit$point.est, fit$se.est,
        ncol(fit$va), ncol(fit$vb), label = labels[i], bc_control = bc_control
      )
    }
  }

  c(list(firth = firth), intervals)
}
