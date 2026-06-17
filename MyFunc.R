logit <- function(prob) {
  log(prob) - log(1 - prob)
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


## Function for checking if two things are equal within numerical precision
same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}


## Functions for wrapping estimation results into a nice format
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
  if (param == "RD") param.est <- linear.predictors
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


## This function is useful for finding limits on the boundary
## It gives 0.5*exp(x)*(-1+sqrt(1+4exp(-x)))
## This is bounded between 0 and 1, and takes value (-1+sqrt(5))/2 at x=0 (Some relation to golden ratio)
## Limits are 0 and 1 as x goes to -infty and +infty respectively
## The function will never return NaN given a numerical input

getPrbAux <- function(x) {
  ifelse((x < 17) & (x > (-500)),
    0.5 * exp(x) * (-1 + (1 + 4 * exp(-x))^0.5),
    ifelse(x < 0, 0, 1)
  )
}
