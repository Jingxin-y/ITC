## 估计用的通用数学、权重变量选择、模型初值和结果辅助函数。
## 原代码：MyFunc.R；generate_RR.R；compare_RR_anchored.R

logit <- function(prob) {
  log(prob) - log(1 - prob)
}

anchored_fit_control <- function(control = list()) {
  defaults <- list(brm_max_step=1000L, brm_optim_maxit=200L, fc_max_step=300L, glm_maxit=200L,
                   brm_thres=1e-6, fc_thres=1e-6, brm_optim_reltol=1e-6)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("Unknown or duplicate fit_control names. Use iteration limits, brm_thres, fc_thres or brm_optim_reltol.")
  iterations <- intersect(names(control),c("brm_max_step","brm_optim_maxit","fc_max_step","glm_maxit"))
  if (any(!vapply(control[iterations],function(x) is.numeric(x) && length(x)==1L && is.finite(x) &&
      x>=1 && x==floor(x) && x<=.Machine$integer.max,logical(1))))
    stop("fit_control iteration limits must be positive integers.")
  tolerances <- intersect(names(control),c("brm_thres","fc_thres","brm_optim_reltol"))
  if (any(!vapply(control[tolerances],function(x) is.numeric(x) && length(x)==1L &&
      is.finite(x) && x>0,logical(1)))) stop("fit_control tolerances must be positive finite numbers.")
  defaults[names(control)] <- control
  defaults[1:4] <- lapply(defaults[1:4],as.integer)
  defaults
}

expit <- function(logodds) {
  1 / (1 + exp(-logodds))
}

log_sum_exp <- function(eta) {
  eta <- as.numeric(eta)
  if (any(!is.finite(eta))) stop("Non-finite linear predictor in log_sum_exp().")

  offset <- max(eta)
  offset + log(sum(exp(eta - offset)))
}

normalize_exp_weights <- function(eta) {
  eta <- as.numeric(eta)
  if (any(!is.finite(eta))) stop("Non-finite linear predictor in weights.")

  eta <- eta - max(eta)
  w <- exp(eta)
  w.mean <- mean(w)
  if (!is.finite(w.mean) || w.mean <= 0) {
    stop("Unable to normalize exponential weights.")
  }

  w / w.mean
}

getlogop <- function(p0, p1) {
  log(p0) + log(p1) - log(1 - p0) - log(1 - p1)
}

getlogrr <- function(p0, p1) {
  log(p1) - log(p0)
}

getatanhrd <- function(p0, p1) {
  atanh(p1 - p0)
}

same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}

WrapResults <- function(point.est, cov, param, name, va, vb, converged, ci.est = NULL) {
  se.est <- sqrt(diag(cov))

  if (is.null(ci.est)) {
    conf.lower <- point.est + stats::qnorm(0.025) * se.est
    conf.upper <- point.est + stats::qnorm(0.975) * se.est
    p.temp <- stats::pnorm(point.est / se.est, 0, 1)
    p.value <- 2 * pmin(p.temp, 1 - p.temp)
  } else {
    conf.lower <- ci.est$low
    conf.upper <- ci.est$up
    p.value <- ci.est$p
  }


  names(point.est) <- names(se.est) <- rownames(cov) <- colnames(cov) <- names(conf.lower) <- names(conf.upper) <- names(p.value) <- name

  coefficients <- cbind(point.est, se.est, conf.lower, conf.upper, p.value)

  linear.predictors <- va %*% point.est[1:ncol(va)]
  if (param == "RR") param.est <- exp(linear.predictors)
  if (param == "RD") param.est <- tanh(linear.predictors)
  if (param == "OR") param.est <- expit(linear.predictors)

  sol <- list(
    param = param, point.est = point.est, se.est = se.est, cov = cov,
    conf.lower = conf.lower, conf.upper = conf.upper, p.value = p.value,
    coefficients = coefficients, param.est = param.est, va = va, vb = vb,
    converged = converged
  )
  class(sol) <- c("brm", "list")
  attr(sol, "hidden") <- c(
    "param", "se.est", "cov", "conf.lower", "conf.upper",
    "p.value", "coefficients", "param.est", "va", "vb", "converged"
  )

  return(sol)
}

getPrbAux <- function(x) {
  ifelse((x < 17) & (x > (-500)),
    0.5 * exp(x) * (-1 + (1 + 4 * exp(-x))^0.5),
    ifelse(x < 0, 0, 1)
  )
}

## 权重协变量与结局模型协变量分开；get.weight() 返回的 gamma 名称可直接传递选择。
weight_covariate_names <- function(data, covariates = NULL, gamma = NULL) {
  if (!is.null(covariates)) return(covariates)
  if (!is.null(names(gamma))) return(names(gamma))
  if (!is.null(gamma) && length(gamma) == 0L) return(character(0))
  covariate_names(data)
}

## 按同一组变量及顺序提取目标均值、方差或协方差。
## 有名称时按名称；无名称时接受原始全维度或已按所选变量排列的输入。
select_weight_summary <- function(x, covariates, all_covariates) {
  if (is.null(x)) return(NULL)
  q <- length(covariates)
  if (is.matrix(x)) {
    if (q == 0L) return(x[integer(0), integer(0), drop = FALSE])
    if (!is.null(rownames(x)) && !is.null(colnames(x))) {
      return(x[covariates, covariates, drop = FALSE])
    }
    idx <- if (nrow(x) == length(all_covariates)) match(covariates, all_covariates) else seq_len(q)
    return(x[idx, idx, drop = FALSE])
  }
  if (q == 0L) return(numeric(0))
  if (!is.null(names(x))) {
    if (!all(covariates %in% names(x))) stop("Target summary is missing selected covariates.")
    return(x[covariates])
  }
  if (length(x) == 1L) return(rep(x, q))
  if (length(x) == length(all_covariates)) return(x[match(covariates, all_covariates)])
  if (length(x) != q) stop("Target summary length must match all or selected covariates.")
  x
}

analysis_terms <- function(data, model_covariates = NULL, effect_modifiers = NULL) {
  if (is.null(model_covariates)) model_covariates <- covariate_names(data)
  interactions <- if (length(effect_modifiers)) paste("x", effect_modifiers, sep = ":") else character(0)
  c("x", model_covariates, interactions)
}

weighted_glm_formula <- function(data, model_covariates = NULL, effect_modifiers = NULL) {
  as.formula(paste("y ~", paste(analysis_terms(data, model_covariates, effect_modifiers), collapse = " + ")))
}

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

make_na_exact <- function(pa = 1L) {
  list(low = rep(NA_real_, pa), up = rep(NA_real_, pa), p = rep(NA_real_, pa))
}

stat_value <- function(x, field) {
  if (is.null(x) || is.null(x[[field]])) return(NA_real_)
  first_value(x[[field]])
}

log_link_start <- function(data, w = NULL, correction = 0.5,
                           max_prob = 0.8, model_covariates = NULL,
                           effect_modifiers = NULL) {
  if (is.null(w)) w <- rep(1, nrow(data))
  terms <- c("(Intercept)", analysis_terms(data, model_covariates, effect_modifiers))
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
  start["(Intercept)"] <- log(p0)
  start[is.na(start)] <- 0
  start
}

add_glm_weights <- function(data, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(data))
  if (length(w) != nrow(data)) {
    stop("Length of w must equal nrow(data).")
  }

  data$.glm_weight <- as.numeric(w)
  data
}

get_glm_stat <- function(fit, stat) {
  if (is.null(fit)) return(NA_real_)

  val <- tryCatch({
    switch(
      stat,
      estimate = unname(stats::coef(fit)["x"]),
      se = unname(summary(fit)$coefficients["x", 2]),
      low = unname(confint.default(fit, level = 0.95)["x", 1]),
      up = unname(confint.default(fit, level = 0.95)["x", 2]),
      p = unname(summary(fit)$coefficients["x", 4]),
      NA_real_
    )
  }, error = function(e) NA_real_)

  first_value(val)
}

exact_safe <- function(param, y, x, va, vb, weight,
                            point.est, se.est, pa, pb, label,
                            max.step = 100,
                            thres = 1e-6,
                            thres.dicho = 1e-3, bc_control = list()) {
  if (any(!is.finite(point.est[seq_len(pa + pb)]))) {
    warning(label, " skipped because point estimates are non-finite.",
            call. = FALSE)
    return(make_na_exact(pa))
  }

  if (any(!is.finite(se.est[seq_len(pa)])) ||
      any(se.est[seq_len(pa)] <= 0)) {
    warning(label, " received non-finite alpha SE; using default exact-search width.",
            call. = FALSE)
    se.est[seq_len(pa)] <- 1
  }

  args <- list(param = param, y = y, x = x, va = va, vb = vb, weight = weight,
               max.step = max.step, thres = thres, thres.dicho = thres.dicho,
               pars = point.est, se = se.est, pa = pa, pb = pb)
  args[names(bc_control)] <- bc_control
  out <- safe_try(do.call(exact, args), label)

  if (is.null(out)) make_na_exact(pa) else out
}


## 原 exact() 通过外部 RcppExports.R 获取此矩阵乘法；此处用等价的 R 运算。
mat_vec_mul <- function(mat, vec) as.vector(mat %*% vec)
