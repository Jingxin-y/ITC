covariate_names <- function(data) {
  grep("^z[0-9]+$", names(data), value = TRUE)
}

covariate_matrix <- function(data) {
  z.cols <- covariate_names(data)
  as.matrix(data[, z.cols, drop = FALSE])
}

analysis_terms <- function(data) {
  c("x", "v.1", covariate_names(data))
}

weighted_glm_formula <- function(data) {
  as.formula(paste("y ~", paste(analysis_terms(data), collapse = " + "), "- 1"))
}

normalize_prob <- function(prob, n_cov, default) {
  if (is.null(prob)) prob <- default
  if (length(prob) == 1) prob <- rep(prob, n_cov)
  as.numeric(prob)
}

binary_covariates <- function(n, prob) {
  z <- sapply(prob, function(p) rbinom(n, 1, p))
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  colnames(z) <- paste0("z", seq_len(ncol(z)))
  z
}

target_prob_vector <- function(target_mean, n_cov) {
  if (is.null(target_mean)) return(NULL)
  if (length(target_mean) == 1) target_mean <- rep(target_mean, n_cov)
  as.numeric(target_mean)
}

rel_ess_binary <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  1 / prod(q^2 / p + (1 - q)^2 / (1 - p))
}

ad_prob_for_ess_theory <- function(ipd_prob,
                                   ess_ratio,
                                   n_cov = length(ipd_prob),
                                   ess_side = "lower",
                                   direction = NULL,
                                   edge_eps = 1e-6) {
  p <- normalize_prob(ipd_prob, n_cov, default = 0.5)
  ess_ratio <- as.numeric(ess_ratio)
  
  if (is.null(direction)) {
    direction <- if (tolower(ess_side) %in% c("upper", "right", "high")) 1 else -1
  }
  if (length(direction) == 1) direction <- rep(direction, length(p))
  direction <- as.numeric(direction)
  
  q_of_s <- function(s) {
    q <- plogis(qlogis(p) + s * direction)
    pmin(pmax(q, edge_eps), 1 - edge_eps)
  }
  
  if (abs(ess_ratio - 1) < 1e-12) {
    return(p)
  }
  
  f <- function(s) rel_ess_binary(p, q_of_s(s)) - ess_ratio
  upper <- 1
  while (f(upper) > 0 && upper < 100) {
    upper <- upper * 2
  }
  
  if (f(upper) > 0) {
    warning("The requested ESS ratio may not be achievable along this direction. Returning the boundary value.",
            call. = FALSE)
    return(q_of_s(upper))
  }
  
  q_of_s(uniroot(f, c(0, upper))$root)
}

effect_vector <- function(intercept, slope, n_cov) {
  c(intercept, rep(slope / sqrt(n_cov), n_cov))
}

scenario_key <- function(event, hypothesis) {
  if (event == "common" && hypothesis == "alternative") {
    "common_alt"
  } else if (event == "common") {
    "common_null"
  } else if (hypothesis == "alternative") {
    "rare_alt"
  } else {
    "rare_null"
  }
}

truth_params <- function(group, event, hypothesis, n_cov) {
  settings <- list(
    IPD = list(
      common_alt = c(alpha = 0.4, beta0 = 1.3, slope = -0.5),
      common_null = c(alpha = 0, beta0 = 1.5, slope = 0.6),
      rare_alt = c(alpha = 0.7, beta0 = -5, slope = -0.5),
      rare_null = c(alpha = 0, beta0 = -4.7, slope = 0.5)
    ),
    AD = list(
      common_alt = c(alpha = 0.7, beta0 = -0.4, slope = -0.5),
      common_null = c(alpha = 0, beta0 = 1.5, slope = 0.6),
      rare_alt = c(alpha = 1.1, beta0 = -4.6, slope = -0.8),
      rare_null = c(alpha = 0, beta0 = -4.7, slope = 0.5)
    )
  )
  
  vals <- settings[[group]][[scenario_key(event, hypothesis)]]
  list(alpha = vals[["alpha"]],
       beta = effect_vector(vals[["beta0"]], vals[["slope"]], n_cov))
}

ad_covariate_prob <- function(target_mean, n_cov,
                              ess_ratio = NULL,
                              ess_side = "lower",
                              ipd_prob = 0.5,
                              ad_prob = NULL,
                              direction = NULL) {
  ess_ratio <- optional_arg(ess_ratio)
  target_mean <- optional_arg(target_mean)
  ad_prob <- optional_arg(ad_prob)
  
  target_mean <- target_prob_vector(target_mean, n_cov)
  if (!is.null(target_mean)) {
    return(target_mean)
  }
  
  if (!is.null(ess_ratio)) {
    return(
      ad_prob_for_ess_theory(
        ipd_prob = ipd_prob,
        ess_ratio = ess_ratio,
        n_cov = n_cov,
        ess_side = ess_side,
        direction = direction
      )
    )
  }
  
  normalize_prob(ad_prob, n_cov, default = 0.35)
}

data.generation <- function(param, group, n, event, hypothesis,
                            ess_ratio = NULL,
                            ess_side = "lower",
                            target_mean = NULL,
                            n_cov = 3,
                            ipd_prob = 0.5,
                            ad_prob = NULL,
                            direction = NULL) {
  n_cov <- as.integer(n_cov)
  target_mean.used <- rep(NA_real_, n_cov)
  
  if (group == "IPD") {
    z.prob <- normalize_prob(ipd_prob, n_cov, default = 0.5)
  } else if (group == "AD") {
    ad.target <- ad_covariate_prob(
      target_mean = target_mean,
      n_cov = n_cov,
      ess_ratio = ess_ratio,
      ess_side = ess_side,
      ipd_prob = ipd_prob,
      ad_prob = ad_prob,
      direction = direction
    )
    z.prob <- ad.target
    target_mean.used <- z.prob
  }
  
  z <- binary_covariates(n, z.prob)
  truth <- truth_params(group, event, hypothesis, n_cov)
  alpha.true <- truth$alpha
  beta.true <- truth$beta
  
  v.1 <- rep(1, n)
  v <- cbind(v.1, z)
  pscore.true <- rep(0.5, n)
  p0p1.true <- getProbRR(as.matrix(v.1, ncol = 1) %*% alpha.true,
                         v %*% beta.true)
  
  x <- rbinom(n, 1, pscore.true)
  pA.true <- p0p1.true[, 1]
  pA.true[x == 1] <- p0p1.true[x == 1, 2]
  y <- rbinom(n, 1, pA.true)
  
  count <- c(
    Na0 = sum(x == 0),
    Na1 = sum(x == 1),
    N0_1 = sum(y[x == 0]),
    N1_1 = sum(y[x == 1])
  )
  
  list(
    data = data.frame(y = y, x = x, v.1 = v.1, z, check.names = FALSE),
    count = count,
    design = list(
      target_mean = target_mean.used,
      ess_ratio = if (is.null(ess_ratio)) NA_real_ else ess_ratio,
      ess_side = ess_side,
      n_cov = n_cov,
      covariate_prob = z.prob
    )
  )
}