## RR: theta=log(RR), phi=log(odds product).
## Derivative source: brmplus commit 3b628a6d32dc7fe0f31fbd3ab181aecd3866f982,
## compare/MLE_Point_Firth_for_RR.R; see BRMPLUS_LICENSE.txt.
## Scalar derivatives reuse the original brmplus RR blocks; va/vb enter only
## in estimation/helpers/firth_augmentation.R. Supports effect modifiers (pa > 1).
firth_rr_derivatives <- function(alpha, beta, x, va, vb) {
  probs <- getProbRR(as.vector(va %*% alpha), as.vector(vb %*% beta))
  p0 <- pmin(pmax(probs[, 1], 1e-8), 1 - 1e-8)
  p1 <- pmin(pmax(probs[, 2], 1e-8), 1 - 1e-8)
  p <- ifelse(x == 1, p1, p0)
  q0 <- 1 - p0
  q1 <- 1 - p1
  denom <- q0 + q1

  ## First and second derivatives of log(p0).
  g.a <- -q0 / denom
  g.b <- q0 * q1 / denom
  g.aa <- ((p0 - p1) * g.a - q0 * p1) / denom^2
  g.ab <- q0 * q1 * (p0 - p1) / denom^3
  g.bb <- -(p0 * q1^2 + p1 * q0^2) * g.b / denom^2

  ## Convert to derivatives of p_x = exp(log(p0) + x*theta).
  list(p = p, d.a = p * (g.a + x), d.b = p * g.b,
       d.aa = p * ((g.a + x)^2 + g.aa),
       d.ab = p * ((g.a + x) * g.b + g.ab),
       d.bb = p * (g.b^2 + g.bb),
       interior = all(is.finite(probs) & probs > 1e-8 & probs < 1 - 1e-8))
}

## RR-only local root rescue, entered only after primary nonconvergence.
## Uses the ORIGINAL FC equation/scale and the REMAINING max.step budget.
## nleqslv is required only on this branch; missing dependency retains failure.
## See estimation/README_rr_fc_rescue_20261019.txt for diagnostics and tests.
.firth_rr_rescue <- function(primary, y, x, va, vb, weight, max.step,
                                    thres, pa, pb, stepmax = 5) {
  if (isTRUE(primary$convergence)) return(primary)
  if (any(!y %in% c(0, 1)) || any(!x %in% c(0, 1)) ||
      any(!is.finite(c(y, x, va, vb, weight))) ||
      length(y) != length(x) || length(weight) != length(y) ||
      nrow(va) != length(y) || nrow(vb) != length(y) ||
      ncol(va) != pa || ncol(vb) != pb ||
      any(weight < 0) || !any(weight > 0)) stop("Invalid RR FC rescue data.")
  if (length(stepmax) != 1L || !is.finite(stepmax) || stepmax <= 0)
    stop("stepmax must be positive and finite.")
  original.diagnostics <- primary$diagnostics
  scale <- original.diagnostics$score_scale
  remaining <- as.integer(max.step - primary$step)
  attempts <- list(); evaluations <- invalid <- iterations <- 0L
  finish.failure <- function(reason) {
    primary$diagnostics$rescue <- list(attempted=length(attempts)>0L,
      accepted=FALSE, reason=reason, attempts=attempts,
      iterations=iterations, evaluations=evaluations, invalid_candidates=invalid,
      maximum_iterations=max(0L,remaining))
    primary$diagnostics$primary_diagnostics <- original.diagnostics
    primary$step <- primary$step + iterations
    primary$diagnostics$iterations <- primary$step
    primary$diagnostics$evaluations <- original.diagnostics$evaluations + evaluations
    primary
  }
  if (remaining < 1L) return(finish.failure("No iterations left in original max.step budget"))
  if (!requireNamespace("nleqslv",quietly=TRUE))
    return(finish.failure("nleqslv is unavailable; original failure retained"))
  if (length(scale) != pa+pb || any(!is.finite(scale)) || any(scale <= 0))
    return(finish.failure("Missing original score scale; original failure retained"))
  w <- weight/max(weight); w <- w/mean(w)
  evaluate <- function(par) {
    if (length(par) != pa+pb || any(!is.finite(par))) stop("Invalid coefficients")
    terms <- firth_rr_derivatives(par[seq_len(pa)],par[pa+seq_len(pb)],x,va,vb)
    if (!isTRUE(terms$interior) || any(!is.finite(terms$p)) ||
        any(terms$p <= 1e-8 | terms$p >= 1-1e-8)) stop("Outside original probability interior")
    z <- firth_score_adjustment(y,va,vb,w,terms)
    if (any(!is.finite(c(z$adjusted.score,z$information,z$covariance))))
      stop("Nonfinite original FC equation or covariance")
    z
  }
  fn <- function(par) {
    evaluations <<- evaluations+1L
    tryCatch(evaluate(par)$adjusted.score/scale,error=function(e) {
      invalid <<- invalid+1L
      rep(NaN,pa+pb)
    })
  }
  # Start from the primary solver's archived best FINITE candidate.
  # Each method uses that same candidate, rather than a failed rescue iterate.
  for (method in c("Newton","Broyden")) {
    allowance <- remaining-iterations
    if (allowance < 1L) break
    fit <- tryCatch(nleqslv::nleqslv(primary$par,fn,method=method,global="none",
      control=list(maxit=allowance,stepmax=stepmax,ftol=min(thres,1e-10),
        xtol=1e-12,allowSingular=TRUE)),error=identity)
    if (inherits(fit,"error")) {
      # Failed calls cannot report how much of their internal iteration budget
      # was consumed. Charge the full allowance and stop to retain the cap.
      iterations <- iterations+allowance
      attempts[[length(attempts)+1L]] <- list(method=method,global="none",
        accepted=FALSE,error=conditionMessage(fit),charged_iterations=allowance)
      break
    }
    iterations <- iterations+fit$iter
    # Recompute from the ORIGINAL binary y and the ORIGINAL score scale.
    # Never accept a termination code, parameter-step test, or residual minimum.
    checked <- tryCatch(evaluate(fit$x),error=function(e)NULL)
    maximum <- if(is.null(checked)) Inf else max(abs(checked$adjusted.score/scale))
    pd <- !is.null(checked) &&
      !inherits(try(chol(checked$information),silent=TRUE),"try-error") &&
      !inherits(try(chol(checked$covariance),silent=TRUE),"try-error")
    accepted <- isTRUE(pd) && is.finite(maximum) && maximum <= thres &&
      iterations <= remaining
    attempts[[length(attempts)+1L]] <- list(method=method,global="none",
      stepmax=stepmax,termcd=fit$termcd,message=fit$message,iterations=fit$iter,
      package_function_evaluations=fit$nfcnt,original_scaled_score=maximum,
      positive_definite=pd,accepted=accepted)
    if (accepted) {
      ans <- primary
      ans$par <- fit$x; ans$convergence <- TRUE
      ans$step <- primary$step+iterations
      ans$diagnostics$primary_diagnostics <- original.diagnostics
      ans$diagnostics$solver <- "joint_damped_newton_then_nleqslv_local_rr"
      ans$diagnostics$reason <- "original adjusted-score tolerance reached by local RR rescue"
      ans$diagnostics$iterations <- ans$step
      ans$diagnostics$evaluations <- original.diagnostics$evaluations+evaluations
      ans$diagnostics$adjusted_score <- checked$adjusted.score
      ans$diagnostics$max_abs_score <- max(abs(checked$adjusted.score))
      ans$diagnostics$max_scaled_score <- maximum
      ans$diagnostics$working_score_scale <- scale
      ans$diagnostics$rescue <- list(attempted=TRUE,accepted=TRUE,
        method=method,global="none",stepmax=stepmax,
        package_version=as.character(utils::packageVersion("nleqslv")),
        iterations=iterations,evaluations=evaluations,invalid_candidates=invalid,
        maximum_iterations=remaining,attempts=attempts)
      return(ans)
    }
  }
  finish.failure("No acceptable original-equation root found within remaining budget")
}

max.likelihood.firth.rr <- function(param, y, x, va, vb, alpha.start, beta.start,
                                    weight, max.step, thres, pa, pb) {
  primary <- max.likelihood.firth(y, x, va, vb, alpha.start, beta.start, weight,
                                  max.step, thres, pa, pb, firth_rr_derivatives)
  # Keep every successful primary return exactly as before, including metadata.
  if (isTRUE(primary$convergence)) return(primary)
  .firth_rr_rescue(primary, y, x, va, vb, weight, max.step, thres, pa, pb)
}
