## Two-arm RR + OP estimation model. logop is not the three-arm log GOP.
## Data generation uses data/helpers/getProbRR_GOP.R independently.

## Exact RR/OP probability mapping, shared by scalar and vector entry points.
## Let a=abs(logrr), r=exp(-a), H=max(p0,p1), L=r*H. Rationalizing
## the quadratic gives H=2/{1+r+sqrt((1-r)^2+4*exp(-a-logop))}.
## Evaluate its positive sums in log space: no subtractive cancellation,
## overflowing positive exponentials, or artificial +/-12 boundary switch.
.rr_op_probabilities <- function(logrr, logop) {
  if (!is.numeric(logrr) || !is.numeric(logop) ||
      is.complex(logrr) || is.complex(logop) ||
      any(!is.finite(logrr)) || any(!is.finite(logop)))
    stop("RR/OP predictors must be finite real numeric values.", call. = FALSE)
  logrr <- as.numeric(logrr); logop <- as.numeric(logop)
  if (!length(logrr) && !length(logop))
    return(matrix(numeric(), 0L, 2L, dimnames = list(NULL, c("p0", "p1"))))
  n <- max(length(logrr), length(logop))
  if (!length(logrr) || !length(logop) ||
      !all(c(length(logrr), length(logop)) %in% c(1L, n)))
    stop("RR/OP predictors must have equal lengths or be scalar.", call. = FALSE)
  logrr <- rep_len(logrr, n); logop <- rep_len(logop, n)
  logadd <- function(u, v) {
    m <- pmax(u, v)
    ans <- m + log1p(exp(pmin(u, v) - m))
    ans[is.infinite(m)] <- m[is.infinite(m)]
    ans
  }
  a <- abs(logrr)
  log.square <- 2 * log(-expm1(-a))
  log.root <- .5 * logadd(log.square, log(4) - a - logop)
  log.denominator <- logadd(log1p(exp(-a)), log.root)
  ## The exact H is <=1; this guard only removes floating-point overshoot.
  log.high <- pmin(0, log(2) - log.denominator)
  high <- exp(log.high)
  low <- exp(log.high - a)
  ## At extremely large finite predictors, machine rounding/underflow may
  ## still produce 1/0. Callers retain their own interior-probability checks.
  cbind(p0 = ifelse(logrr >= 0, low, high),
        p1 = ifelse(logrr >= 0, high, low))
}

#' Calculate two-arm risks from log RR and log odds product
#'
#' @param logrr Log relative risk, or a length-two c(logrr, logop) vector.
#' @param logop Log odds product.
#' @details OP=p0*p1/((1-p0)*(1-p1)); RR=p1/p0. Uses an exact,
#' numerically stable mapping without a parameter boundary approximation.
#' @return A length-two vector c(p0, p1).
#' @examples getProbScalarRR(0, 0)
#' @export
getProbScalarRR <- function(logrr, logop = NA) {
  if (length(logrr) == 2L) {
    logop <- logrr[2L]; logrr <- logrr[1L]
  }
  if (length(logrr) != 1L || length(logop) != 1L)
    stop("getProbScalarRR requires one RR/OP pair.", call. = FALSE)
  as.vector(.rr_op_probabilities(logrr, logop))
}

#' Calculate two-arm risks from log RR and log odds product (vectorized)
#'
#' @param logrr Log relative risks; also accepts a two-column predictor matrix
#' or a length-two c(logrr, logop) vector when logop is omitted.
#' @param logop Log odds products. A scalar is recycled to the common length.
#' @return A matrix with columns p0 and p1, including for a single pair.
#' @examples getProbRR(c(0, log(2)), c(-13, 0))
#' @export
getProbRR <- function(logrr, logop = NA) {
  if (is.matrix(logrr) && ncol(logrr) == 2L) {
    logop <- logrr[, 2L]; logrr <- logrr[, 1L]
  } else if (length(logop) == 1L && is.na(logop) && length(logrr) == 2L) {
    logop <- logrr[2L]; logrr <- logrr[1L]
  }
  .rr_op_probabilities(logrr, logop)
}
