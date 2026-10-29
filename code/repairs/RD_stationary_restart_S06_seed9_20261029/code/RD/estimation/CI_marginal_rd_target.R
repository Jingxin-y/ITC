## RD-specific adaptation of the existing RR Blaker-type bootstrap inversion.
## Both the candidate and FIXED observed target are bootstrapped, including FC anchors.
## Conditional on covariates, treatment, fitting/standardization weights; no exact
## finite-sample coverage guarantee. The original RR source is unchanged.
## Saturated effect and OP designs are supported. All effect outputs are natural RD.

.rd_exact_risks <- function(theta, phi) {
  z <- rd_probabilities(theta, phi)
  s0 <- z$p0 * (1 - z$p0); s1 <- z$p1 * (1 - z$p1)
  den <- pmax(s0 + s1, 1e-300); q <- 1 - tanh(theta)^2
  da0 <- -s0 / den * q; da1 <- s1 / den * q; db <- s0 * s1 / den
  list(p0 = z$p0, p1 = z$p1, logp0 = log(z$p0), logp1 = log(z$p1),
    t0 = da0 / pmax(z$p0, 1e-300), t1 = da1 / pmax(z$p1, 1e-300),
    s0 = db / pmax(z$p0, 1e-300), s1 = db / pmax(z$p1, 1e-300))
}

.rd_exact_weights <- function(w, n, label) {
  if (is.null(w)) w <- rep(1, n)
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w)) ||
      any(w < 0) || !any(w > 0))
    stop(label, " must be finite nonnegative weights with positive sum.", call. = FALSE)
  w <- w / max(w)
  w / sum(w)
}

## Deferred boundary-rescue arithmetic. Keep the regular successful path intact.
## Compute small risks AND their complements directly; never subtract a risk
## rounded near one, and keep 1-|tanh(theta)| on its logistic tail scale.
.rd_boundary_logsum <- function(a, b) {
  top <- pmax(a, b)
  ans <- top + log1p(exp(-abs(a-b)))
  ans[is.infinite(top) & top < 0] <- -Inf
  ans
}
.rd_boundary_log1mexp <- function(a) {
  out <- numeric(length(a)); small <- a < -log(2)
  out[small] <- log1p(-exp(a[small]))
  out[!small] <- log(-expm1(a[!small]))
  out
}
.rd_boundary_risks <- function(theta, phi) {
  d <- tanh(theta); t <- abs(d)
  log.gap <- log(2) + plogis(-2*abs(theta), log.p=TRUE)
  log.t <- log(t)
  log.minor <- function(op) {
    u <- exp(-abs(op))
    a <- t * (-expm1(-abs(op))); b <- 2*exp(-abs(op)/2)
    scale <- pmax(a,b); root <- numeric(length(theta)); active <- scale > 0
    root[active] <- scale[active]*sqrt((a[active]/scale[active])^2 + (b[active]/scale[active])^2)
    out <- numeric(length(theta)); zero <- t == 0
    out[zero] <- plogis(op[zero]/2, log.p=TRUE)
    high <- op >= 0 & !zero
    out[high] <- log(2) + log.gap[high] -
      log(2-t[high]+t[high]*u[high]+root[high])
    low <- op < 0 & !zero
    out[low] <- log(2) + op[low] + log.gap[low] -
      log(u[low]*(2-t[low])+t[low]+root[low])
    out
  }
  lm <- log.minor(phi); lc <- log.minor(-phi)
  # Complement-and-arm reversal changes phi to -phi while preserving RD.
  positive <- theta >= 0
  lp0 <- ifelse(positive, lm, .rd_boundary_logsum(log.t,lm))
  lp1 <- ifelse(positive, .rd_boundary_logsum(log.t,lm), lm)
  lq0 <- ifelse(positive, .rd_boundary_logsum(log.t,lc), lc)
  lq1 <- ifelse(positive, lc, .rd_boundary_logsum(log.t,lc))
  normalize.pair <- function(lp,lq) {
    # The smaller probability is evaluated directly; obtain its large complement
    # with log1mexp, avoiding cancellation on either side of RD=0.
    p.small <- lp <= lq
    lp[!p.small] <- .rd_boundary_log1mexp(lq[!p.small])
    lq[p.small] <- .rd_boundary_log1mexp(lp[p.small])
    list(lp=lp,lq=lq)
  }
  p0 <- normalize.pair(lp0,lq0); p1 <- normalize.pair(lp1,lq1)
  gap <- exp(log.gap)
  list(p0=exp(p0$lp),p1=exp(p1$lp),q0=exp(p0$lq),q1=exp(p1$lq),
       logp0=p0$lp,logp1=p1$lp,logq0=p0$lq,logq1=p1$lq,
       theta.derivative=gap*(2-gap))
}

## Internal constrained optimizer; returned closures also allow independent
## numerical verification without running a bootstrap interval.
.rd_exact_problem <- function(y, x, va, vb, weight, pars,
    va.target, vb.target, standardization.weight, intercept,
    alpha.bound, beta.bound, optim.maxit, optim.reltol, score.tol,
    root.tol, root.maxit, constraint.tol, eps, check.deadline,
    optim.retries = 2L, optim.retry.maxit = max(1000L, optim.maxit),
    optim.multistart = TRUE, optim.multistart.max = 3L, optim.fallback = TRUE,
    optim.boundary.rescue = TRUE) {
  pa <- ncol(va); pb <- ncol(vb); rest <- setdiff(seq_len(pa), intercept)
  fw <- .rd_exact_weights(weight, length(y), "weight") * length(y)
  omega <- .rd_exact_weights(standardization.weight, nrow(va.target),
                                   "standardization.weight")
  start <- c(pars[rest], pars[pa + seq_len(pb)])
  bound <- c(rep(alpha.bound, length(rest)), rep(beta.bound, pb))
  ## Repair only the optimizer start; pars still defines the observed target.
  outside <- abs(start) > bound
  start[outside] <- sign(start[outside]) *
    (bound[outside] - pmin(1e-6, bound[outside]/2))
  target <- function(p) marginal_rd_target(p, va.target, vb.target, omega, TRUE)
  stable.boundary <- FALSE
  stable.target.gradient <- function(par) {
    eta <- drop(va.target %*% par[seq_len(pa)])
    gap <- 2*plogis(-2*abs(eta))
    c(colSums(va.target * (omega * gap*(2-gap))), rep(0,pb))
  }
  unpack <- function(nuisance, psi) {
    check.deadline("marginal constraint")
    alpha <- numeric(pa); alpha[rest] <- nuisance[seq_along(rest)]
    beta <- nuisance[length(rest) + seq_len(pb)]
    offset <- drop(va.target %*% alpha)
    active <- omega > 0
    ## atanh(mean RD) lies between the smallest/largest conditional atanh RD.
    ## This supplies an exact root bracket, independent of beta; no intercept cap.
    lower <- atanh(psi) - max(offset[active]); upper <- atanh(psi) - min(offset[active])
    f <- function(a0) {
      alpha[intercept] <- a0
      marginal_rd_target(c(alpha, beta), va.target, vb.target, omega) - psi
    }
    if (stable.boundary) {
      # Solve the same constraint on a relative tail-probability scale.
      # Absolute RD residuals near +/-1 lose significant digits and used to
      # accept a visibly inaccurate intercept before likelihood optimization.
      sign.psi <- if (psi >= 0) 1 else -1
      log.gap <- log1p(-abs(psi))
      f <- function(a0) {
        eta <- a0 + offset
        pieces <- log(omega[active]) + log(2) +
          plogis(-2*sign.psi*eta[active], log.p=TRUE)
        top <- max(pieces)
        log.mean.gap <- top + log(sum(exp(pieces-top)))
        -sign.psi*(log.mean.gap-log.gap)
      }
    }
    if (upper - lower <= root.tol) {
      a0 <- (lower + upper)/2
    } else {
      fl <- f(lower); fu <- f(upper)
      if (abs(fl) <= root.tol) a0 <- lower else if (abs(fu) <= root.tol) a0 <- upper else {
        if (!is.finite(fl) || !is.finite(fu) || fl > 0 || fu < 0)
          stop("Could not bracket the marginal intercept constraint.")
        a0 <- uniroot(f, c(lower, upper), f.lower = fl, f.upper = fu,
                      tol = root.tol, maxiter = root.maxit, check.conv = TRUE)$root
      }
    }
    alpha[intercept] <- a0
    c(alpha, beta)
  }
  evaluate <- function(nuisance, psi, y.local) {
    par <- unpack(nuisance, psi)
    if (stable.boundary) {
      z <- .rd_boundary_risks(drop(va %*% par[seq_len(pa)]),
                              drop(vb %*% par[pa+seq_len(pb)]))
      lp <- ifelse(x==1,z$logp1,z$logp0); lq <- ifelse(x==1,z$logq1,z$logq0)
      clipped <- lp <= log(eps) | lq <= log(eps)
      # Retain the original eps-clipped likelihood and its zero derivative
      # outside the clipping region, using accurate logs to evaluate it.
      llp <- pmin(log1p(-eps),pmax(log(eps),lp))
      llq <- pmin(log1p(-eps),pmax(log(eps),lq))
      nll <- -sum(fw*ifelse(y.local==1,llp,llq))
      s0 <- z$p0*z$q0; s1 <- z$p1*z$q1; den <- pmax(s0+s1,1e-300)
      residual <- ifelse(y.local==1,exp(lq),-exp(lp))
      sa <- residual*ifelse(x==1,1,-1)*z$theta.derivative/den
      sb <- residual*ifelse(x==1,s0,s1)/den
      sa[clipped] <- sb[clipped] <- 0
      full.grad <- -c(colSums(va*(fw*sa)),colSums(vb*(fw*sb)))
      tg <- target(par); tg$gradient <- stable.target.gradient(par)
      free <- c(rest,pa+seq_len(pb))
      if (!is.finite(tg$gradient[intercept]) || tg$gradient[intercept] <= 0)
        stop("Nonpositive marginal-target derivative for the intercept.")
      grad <- full.grad[free] - full.grad[intercept]*tg$gradient[free]/tg$gradient[intercept]
      return(list(value=nll,gradient=grad,par=par,target=tg$value,
        constraint.error=tg$value-psi,probabilities=exp(lp),clipped.probabilities=sum(clipped)))
    }
    z <- .rd_exact_risks(drop(va %*% par[seq_len(pa)]),
                              drop(vb %*% par[pa + seq_len(pb)]))
    raw <- ifelse(x == 1, z$p1, z$p0)
    p <- pmin(pmax(raw, eps), 1-eps)
    nll <- -sum(fw * (y.local*log(p) + (1-y.local)*log1p(-p)))
    scale <- fw * (y.local-p)/(1-p)
    scale[raw <= eps | raw >= 1-eps] <- 0
    full.grad <- -c(colSums(va * (scale*ifelse(x == 1, z$t1, z$t0))),
                    colSums(vb * (scale*ifelse(x == 1, z$s1, z$s0))))
    tg <- target(par)
    free <- c(rest, pa + seq_len(pb))
    if (!is.finite(tg$gradient[intercept]) || tg$gradient[intercept] <= 0)
      stop("Nonpositive marginal-target derivative for the intercept.")
    grad <- full.grad[free] - full.grad[intercept]*tg$gradient[free]/tg$gradient[intercept]
    list(value = nll, gradient = grad, par = par, target = tg$value,
         constraint.error = tg$value-psi, probabilities = p,
         clipped.probabilities = sum(raw <= eps | raw >= 1-eps))
  }
  ## Bounded warm-start storage: never keep a history of bootstrap responses.
  ## It changes optimizer starts only, not the response, likelihood or constraint.
  observed.starts <- list(); previous <- NULL
  cell.means.start <- function(yy) {
    key <- do.call(paste, c(as.data.frame(cbind(va, vb)), sep = ":"))
    rows <- which(!duplicated(key)); cell <- match(key, key[rows])
    ca <- va[rows, , drop = FALSE]; cb <- vb[rows, , drop = FALSE]
    qa <- qr(ca); qb <- qr(cb)
    if (qa$rank < pa || qb$rank < pb) return(NULL)
    risk <- matrix(NA_real_, length(rows), 2L)
    for (j in seq_along(rows)) for (arm in 0:1) {
      take <- cell == j & x == arm; total <- sum(fw[take])
      if (total <= 0) return(NULL)
      ## Half-event smoothing is ONLY an optimizer start. All objective/score
      ## evaluations continue to use the original yy and unchanged weights.
      risk[j, arm+1L] <- (sum(fw[take]*yy[take]) + .5)/(total + 1)
    }
    a <- qr.coef(qa, atanh(risk[,2L]-risk[,1L]))
    b <- qr.coef(qb, qlogis(risk[,1L])+qlogis(risk[,2L]))
    z <- c(a[rest], b)
    if (any(!is.finite(z))) NULL else z
  }
  profile <- function(psi, y.local = y) {
    check.deadline("profile fit")
    stable.boundary <<- FALSE
    last <- NULL; last.par <- NULL
    get <- function(v) {
      check.deadline("profile objective")
      if (is.null(last.par) || !identical(v, last.par)) {
        last <<- evaluate(v, psi, y.local); last.par <<- v
      }
      last
    }
    ## Retry the SAME response and constraint. No RNG calls, resampling or
    ## relaxed score/constraint thresholds are used to rescue a numerical fit.
    next.start <- start; best <- NULL; attempts <- list(); tried <- list()
    retry.attempts <- multistart.attempts <- fallback.attempts <- boundary.attempts <- 0L
    run.optimizer <- function(initial, method, label, limit, reltol, pgtol) {
      check.deadline(paste("profile", method, label))
      initial <- pmin(pmax(unname(initial), -bound), bound)
      tried[[length(tried)+1L]] <<- initial
      last <<- NULL; last.par <<- NULL
      raw <- NULL
      t0 <- as.numeric(Sys.time())
      fit <- tryCatch({
        fn <- function(v) get(v)$value/length(y)
        gr <- function(v) get(v)$gradient/length(y)
        if (method == "L-BFGS-B") {
          raw <- optim(initial, fn, gr, method = method, lower = -bound, upper = bound,
            control = list(maxit = limit, factr = reltol/.Machine$double.eps, pgtol = pgtol))
        } else {
          raw <- nlminb(initial, fn, gradient = gr, lower = -bound, upper = bound,
            control = list(iter.max = limit, eval.max = max(2000L, 5*limit),
                           rel.tol = max(1e-14, reltol), x.tol = 1e-10))
          raw$counts <- raw$evaluations
        }
        z <- get(raw$par)
        pg <- z$gradient/length(y)
        pg[(raw$par <= -bound+1e-7 & pg > 0) | (raw$par >= bound-1e-7 & pg < 0)] <- 0
        ok <- raw$convergence == 0L && all(is.finite(c(z$value,z$par,pg))) &&
          abs(z$constraint.error) <= constraint.tol && max(abs(pg)) <= score.tol
        list(converged = ok, psi = psi, par = z$par, alpha = z$par[seq_len(pa)],
             beta = z$par[pa+seq_len(pb)], value = z$value, target = z$target,
             constraint.error = z$constraint.error, optimizer.code = raw$convergence,
             max.projected.score = max(abs(pg)), counts = raw$counts,
             nuisance.boundary = any(abs(raw$par) >= bound-1e-6),
             clipped.probabilities = z$clipped.probabilities,
              optimizer.method = method, start.label = label,
              arithmetic = if(stable.boundary) "stable_boundary" else "original",
              optimizer.message = if(is.null(raw$message)) "" else raw$message,
             reason = if (ok) "converged" else "optimizer, projected score or constraint tolerance not satisfied")
      }, error = function(e) {
        if (inherits(e, "marginal_exact_timeout")) stop(e)
        list(converged = FALSE, psi = psi, reason = conditionMessage(e),
             par = if (is.null(last)) pars else last$par, value = NA_real_,
              constraint.error = NA_real_, optimizer.code = 99L,
              optimizer.method = method, start.label = label)
      })
      attempts[[length(attempts)+1L]] <<- list(attempt=length(attempts)+1L,
        method=method, start.label=label, start=initial, maxit=limit, reltol=reltol,pgtol=pgtol,
        converged=fit$converged, optimizer.code=fit$optimizer.code,
        arithmetic=if(stable.boundary) "stable_boundary" else "original",
        optimizer.message=fit$optimizer.message,
        value=fit$value, max.projected.score=fit$max.projected.score,
        constraint.error=fit$constraint.error, reason=fit$reason,
        elapsed.sec=as.numeric(Sys.time())-t0)
      ## A failed low objective must never replace a verified converged fit.
      if (is.null(best) || (isTRUE(fit$converged) && !isTRUE(best$converged)) ||
          (identical(isTRUE(fit$converged), isTRUE(best$converged)) &&
           is.finite(fit$value) && (!is.finite(best$value) || fit$value < best$value))) best <<- fit
      ## Continue from the returned nuisance coefficients, or the last valid
      ## objective evaluation after an error. Bounds remain fixed.
      resume <- if (!is.null(raw) && all(is.finite(raw$par))) raw$par else last.par
      if (length(resume) == length(start) && all(is.finite(resume)))
        next.start <<- pmin(pmax(resume, -bound), bound)
      fit
    }
    for (attempt in seq_len(optim.retries + 1L)) {
      limit <- if (attempt == 1L) optim.maxit else optim.retry.maxit
      reltol <- max(min(.Machine$double.eps, optim.reltol), optim.reltol/100^(attempt-1L))
      pgtol <- score.tol/10/100^(attempt-1L)
      z <- run.optimizer(next.start, "L-BFGS-B", if (attempt == 1L) "original" else "resume",
                         limit, reltol, pgtol)
      retry.attempts <- retry.attempts+1L
      if (isTRUE(z$converged)) break
    }
    alternatives <- list()
    add.start <- function(v, label) {
      if (is.null(v) || length(v) != length(start) || any(!is.finite(v))) return(invisible(NULL))
      v <- pmin(pmax(unname(v), -bound), bound)
      all.starts <- c(tried, lapply(alternatives, function(z) z$par))
      if (!any(vapply(all.starts, function(z) max(abs(v-z)) < 1e-10, logical(1))))
        alternatives[[length(alternatives)+1L]] <<- list(par=v, label=label)
      invisible(NULL)
    }
    if (!isTRUE(best$converged) && optim.multistart) {
      if (!is.null(previous) && identical(previous$response, y.local))
        add.start(previous$nuisance, "same_response_neighbor")
      if (length(observed.starts)) {
        nearest <- which.min(vapply(observed.starts, function(z) abs(z$psi-psi), numeric(1)))
        add.start(observed.starts[[nearest]]$nuisance, "observed_neighbor")
      }
      check.deadline("cell-means optimizer start")
      add.start(cell.means.start(y.local), "cell_means")
      add.start(rep(0, length(start)), "neutral")
      if (length(alternatives) > optim.multistart.max)
        alternatives <- alternatives[seq_len(optim.multistart.max)]
      ## Compare all verified solutions obtained from the additional starts.
      for (initial in alternatives) {
        run.optimizer(initial$par, "L-BFGS-B", initial$label, optim.retry.maxit, reltol, pgtol)
        multistart.attempts <- multistart.attempts+1L
      }
    }
    if (!isTRUE(best$converged) && optim.fallback) {
      resume <- if (!is.null(best$par) && all(is.finite(best$par)))
        best$par[c(rest, pa+seq_len(pb))] else next.start
      fallback.starts <- c(list(list(par=resume, label="best_resume")), alternatives)
      for (initial in fallback.starts) {
        z <- run.optimizer(initial$par, "nlminb", initial$label, optim.retry.maxit, reltol, pgtol)
        fallback.attempts <- fallback.attempts+1L
        if (isTRUE(z$converged)) break
      }
    }
    ## Retry the same stationary candidate with the initial optimizer controls.
    ## Keep the original acceptance tests, likelihood and parameter bounds.
    if (!isTRUE(best$converged) &&
        all(is.finite(c(best$par, best$value, best$max.projected.score, best$constraint.error))) &&
        length(best$par) == pa + pb &&
        length(best$max.projected.score) == 1L && length(best$constraint.error) == 1L &&
        best$max.projected.score <= score.tol && abs(best$constraint.error) <= constraint.tol) {
      run.optimizer(best$par[c(rest, pa + seq_len(pb))], "L-BFGS-B",
                    "stationary_solution_restart", optim.retry.maxit,
                    optim.reltol, score.tol/10)
    }
    if (!isTRUE(best$converged) && optim.boundary.rescue && abs(psi) >= .99) {
      check.deadline("stable boundary rescue")
      stable.boundary <<- TRUE
      # No new response, RNG use, altered bounds or relaxed score/constraint
      # thresholds. Only the arithmetic for the same constrained problem changes.
      resume <- if (!is.null(best$par) && all(is.finite(best$par)))
        best$par[c(rest,pa+seq_len(pb))] else next.start
      rescue.starts <- c(list(list(par=resume,label="stable_boundary_resume")), alternatives,
                         list(list(par=start,label="stable_boundary_original")))
      for (initial in rescue.starts) {
        z <- run.optimizer(initial$par,"L-BFGS-B",initial$label,optim.retry.maxit,
                           optim.reltol,score.tol/10)
        boundary.attempts <- boundary.attempts+1L
        if (isTRUE(z$converged)) break
        if (optim.fallback) {
          z <- run.optimizer(next.start,"nlminb",initial$label,optim.retry.maxit,
                             optim.reltol,score.tol/10)
          boundary.attempts <- boundary.attempts+1L
          if (isTRUE(z$converged)) break
        }
      }
    }
    if (isTRUE(best$converged)) {
      nuisance <- best$par[c(rest, pa+seq_len(pb))]
      previous <<- list(response=y.local, psi=psi, nuisance=nuisance)
      if (identical(y.local, y)) {
        observed.starts[[length(observed.starts)+1L]] <<- list(psi=psi, nuisance=nuisance)
        if (length(observed.starts) > 16L) observed.starts <<- tail(observed.starts, 16L)
      }
    }
    best$optimizer.attempts <- length(attempts)
    best$retry.attempts <- retry.attempts
    best$multistart.attempts <- multistart.attempts
    best$fallback.attempts <- fallback.attempts
    best$boundary.attempts <- boundary.attempts
    best$boundary.recovered <- boundary.attempts > 0L && isTRUE(best$converged)
    best$retry.recovered <- length(attempts) > 1L && isTRUE(best$converged)
    if (length(attempts) > 1L) best$attempts <- attempts
    best
  }
  list(profile = profile, evaluate = evaluate, unpack = unpack, target = target,
       start = start, bound = bound, fit.weights = fw)
}

## Final runtime: RD_final_B600_F1000_20261024; corrected cutoff and B=600/1000/1000.
exact_marginal_rd_target <- function(param = "RD", y, x, va, vb, weight = NULL,
    max.step = 200L, pars, pa = ncol(va), pb = ncol(vb),
    va.target = va, vb.target = vb, standardization.weight = NULL,
    vcov.pars = NULL, psi.se = NULL, level = .95, psi.null = 0,
    optim.maxit = max.step, optim.reltol = 1e-10, score.tol = 1e-4,
    thres.dicho = 1e-3, bisection.max.step = 15L, grid.mult = 2,
    nsim.coarse = 600L, nsim.confirm = 1000L, nsim.refine = 1000L,
    nsim.p.value = 500L, search.mult = NULL, psi.range = NULL,
    alpha.bound = 15, beta.bound = 15, intercept = 1L,
    root.tol = 1e-10, root.maxit = 200L, constraint.tol = 1e-7,
    time.limit.sec = 10800, seed = NULL, compute.p.value = TRUE,
    compute.interval = TRUE, psi.grid = numeric(0),
    optim.retries = 2L, optim.retry.maxit = max(1000L, optim.maxit),
    optim.multistart = TRUE, optim.multistart.max = 3L, optim.fallback = TRUE,
    optim.boundary.rescue = TRUE) {
  if (!identical(param, "RD")) stop("Only natural marginal RD is supported.", call. = FALSE)
  va <- as.matrix(va); vb <- as.matrix(vb)
  va.target <- as.matrix(va.target); vb.target <- as.matrix(vb.target)
  n <- length(y)
  if (!is.numeric(y) || !is.numeric(x) || n < 2 || length(x) != n ||
      any(!is.finite(c(y,x))) || any(!y %in% 0:1) || any(!x %in% 0:1) ||
      length(unique(x)) != 2L) stop("y/x must be finite binary vectors with both arms.")
  if (!is.numeric(va) || !is.numeric(vb) || nrow(va) != n || nrow(vb) != n ||
      pa != ncol(va) || pb != ncol(vb) || pa < 1 || pb < 1 ||
      any(!is.finite(c(va,vb)))) stop("Invalid fitting designs or pa/pb.")
  if (!is.numeric(va.target) || !is.numeric(vb.target) || nrow(va.target) < 1L ||
      nrow(va.target) != nrow(vb.target) || ncol(va.target) != pa || ncol(vb.target) != pb ||
      any(!is.finite(c(va.target,vb.target)))) stop("Invalid standardization designs.")
  scalar <- function(v, name, positive = TRUE, integer = FALSE) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        (positive && v <= 0) || (integer && v != floor(v))) stop("Invalid ",name,".")
  }
  scalar(optim.retries,"optim.retries",positive=FALSE,integer=TRUE)
  if (optim.retries < 0 || optim.retries > 10) stop("optim.retries must be between 0 and 10.")
  scalar(optim.multistart.max,"optim.multistart.max",integer=TRUE)
  if (optim.multistart.max > 8L) stop("optim.multistart.max must be between 1 and 8.")
  for (name in c("optim.multistart","optim.fallback","optim.boundary.rescue"))
    if (!is.logical(get(name)) || length(get(name)) != 1L || is.na(get(name))) stop("Invalid ",name,".")
  for (name in c("optim.maxit","optim.retry.maxit","max.step","bisection.max.step","root.maxit",
                 "nsim.coarse","nsim.confirm","nsim.refine","nsim.p.value","intercept"))
    scalar(get(name),name,integer=TRUE)
  for (name in c("optim.reltol","score.tol","thres.dicho","root.tol",
                 "constraint.tol","alpha.bound","beta.bound","level","grid.mult"))
    scalar(get(name),name)
  if (level >= 1 || grid.mult <= 1 || intercept > pa) stop("Invalid level/grid.mult/intercept.")
  if (any(abs(va[,intercept]-1) > 1e-12) || any(abs(va.target[,intercept]-1) > 1e-12))
    stop("The specified alpha intercept column must be one in both designs.")
  if (!is.numeric(pars) || length(pars) != pa+pb || any(!is.finite(pars))) stop("Invalid pars.")
  if (!is.numeric(psi.grid) || any(!is.finite(psi.grid))) stop("psi.grid must be finite.")
  scalar(psi.null,"psi.null",positive=FALSE)
  if (abs(psi.null) >= 1 || any(abs(psi.grid) >= 1)) stop("RD null/grid must be strictly in (-1,1).")
  if (!is.numeric(time.limit.sec) || length(time.limit.sec) != 1L ||
      is.na(time.limit.sec) || time.limit.sec <= 0) stop("Invalid time.limit.sec.")
  for (name in c("compute.interval","compute.p.value"))
    if (!is.logical(get(name)) || length(get(name)) != 1L || is.na(get(name))) stop("Invalid ",name,".")
  fw <- .rd_exact_weights(weight,n,"weight")
  if (any(vapply(0:1,function(a) sum(fw[x==a]) <= 0,logical(1))))
    stop("Both treatment arms must have positive fitting weight.")
  if (is.null(standardization.weight)) {
    if (nrow(va.target) != n) stop("Supply standardization.weight for a separate target sample.")
    standardization.weight <- fw
  }
  omega <- .rd_exact_weights(standardization.weight,nrow(va.target),"standardization.weight")
  if (!is.null(seed)) {
    scalar(seed,"seed",positive=FALSE,integer=TRUE)
    if (seed < 0 || seed > .Machine$integer.max) stop("Invalid seed.")
    had.seed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
    if (had.seed) old.seed <- get(".Random.seed",.GlobalEnv)
    on.exit(if (had.seed) assign(".Random.seed",old.seed,.GlobalEnv) else
      if (exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
    set.seed(seed)
  }
  started <- as.numeric(Sys.time())
  check.deadline <- function(stage) {
    if (as.numeric(Sys.time())-started >= time.limit.sec)
      stop(structure(list(message=paste("exact_marginal_rd_target() time limit exceeded at",stage),
                          call=NULL,stage=stage),class=c("marginal_exact_timeout","error","condition")))
  }
  problem <- .rd_exact_problem(y,x,va,vb,fw,pars,va.target,vb.target,omega,
    intercept,alpha.bound,beta.bound,optim.maxit,optim.reltol,score.tol,
    root.tol,root.maxit,constraint.tol,1e-12,check.deadline,
    optim.retries=optim.retries,optim.retry.maxit=optim.retry.maxit,
    optim.multistart=optim.multistart,optim.multistart.max=optim.multistart.max,
    optim.fallback=optim.fallback,optim.boundary.rescue=optim.boundary.rescue)
  target <- problem$target(pars); psi.hat <- target$value
  if (!is.finite(psi.hat)) stop("Nonfinite marginal RD at pars.")
  psi.se.source <- "not used"
  if (!is.null(psi.se)) {
    scalar(psi.se,"psi.se"); psi.se.source <- "user supplied, search range only"
  } else if (!is.null(vcov.pars)) {
    V <- as.matrix(vcov.pars)
    if (!is.numeric(V) || !all(dim(V)==c(pa+pb,pa+pb)) || any(!is.finite(V)) ||
        max(abs(V-t(V))) > 1e-8*max(1,max(abs(V)))) stop("vcov.pars must be a symmetric full parameter covariance.")
    variance <- drop(crossprod(target$gradient,V%*%target$gradient))
    if (!is.finite(variance) || variance <= 0) stop("Nonpositive delta variance for search range.")
    psi.se <- sqrt(variance); psi.se.source <- "full-parameter delta SE, search range only"
  }
  if (is.null(search.mult)) search.mult <- if(n < 300) 8 else 4
  scalar(search.mult,"search.mult")
  if (is.null(psi.range) && compute.interval) {
    if (is.null(psi.se)) stop("Supply psi.se, full vcov.pars or an explicit psi.range.")
    psi.range <- pmin(1-1e-8, pmax(-1+1e-8, psi.hat + c(-1,1)*search.mult*psi.se))
  }
  if (!is.null(psi.range) && (!is.numeric(psi.range) || length(psi.range) != 2L ||
      any(!is.finite(psi.range)) || any(abs(psi.range) >= 1) || psi.range[1] >= psi.hat || psi.range[2] <= psi.hat))
    stop("psi.range must strictly contain the marginal point estimate.")

  counts <- list(profile.fits=0L,failed.profile.fits=0L,bootstrap.draws=0L,
                 failed.bootstrap.draws=0L,candidate.cache.hits=0L,nuisance.boundary.fits=0L,
                 optimizer.attempts=0L,retried.profile.fits=0L,recovered.profile.fits=0L,
                 multistart.profile.fits=0L,fallback.profile.fits=0L,
                 multistart.attempts=0L,fallback.attempts=0L,
                 recovered.by.multistart=0L,recovered.by.fallback=0L,
                 boundary.rescue.attempts=0L,recovered.by.boundary.rescue=0L)
  failures <- retries <- list(); cache <- new.env(parent=emptyenv()); history <- list()
  fit.profile <- function(psi, yy, stage) {
    z <- problem$profile(psi,yy); counts$profile.fits <<- counts$profile.fits+1L
    attempts <- if (is.null(z$optimizer.attempts)) 1L else z$optimizer.attempts
    counts$optimizer.attempts <<- counts$optimizer.attempts+attempts
    counts$multistart.attempts <<- counts$multistart.attempts+z$multistart.attempts
    counts$fallback.attempts <<- counts$fallback.attempts+z$fallback.attempts
    counts$boundary.rescue.attempts <<- counts$boundary.rescue.attempts+z$boundary.attempts
    if (isTRUE(z$boundary.recovered))
      counts$recovered.by.boundary.rescue <<- counts$recovered.by.boundary.rescue+1L
    if (z$multistart.attempts > 0L) {
      counts$multistart.profile.fits <<- counts$multistart.profile.fits+1L
      if (isTRUE(z$converged) && !isTRUE(z$boundary.recovered) && identical(z$optimizer.method,"L-BFGS-B"))
        counts$recovered.by.multistart <<- counts$recovered.by.multistart+1L
    }
    if (z$fallback.attempts > 0L) {
      counts$fallback.profile.fits <<- counts$fallback.profile.fits+1L
      if (isTRUE(z$converged) && !isTRUE(z$boundary.recovered))
        counts$recovered.by.fallback <<- counts$recovered.by.fallback+1L
    }
    if (attempts > 1L) {
      counts$retried.profile.fits <<- counts$retried.profile.fits+1L
      if (isTRUE(z$converged)) counts$recovered.profile.fits <<- counts$recovered.profile.fits+1L
      if (length(retries) < 25L) retries[[length(retries)+1L]] <<- list(stage=stage,response=yy,fit=z)
    }
    if (isTRUE(z$nuisance.boundary)) counts$nuisance.boundary.fits <<- counts$nuisance.boundary.fits+1L
    if (!isTRUE(z$converged)) {
      counts$failed.profile.fits <<- counts$failed.profile.fits+1L
      if (length(failures) < 25L) failures[[length(failures)+1L]] <<- list(stage=stage,response=yy,fit=z)
    }
    z
  }
  candidate <- function(psi) {
    key <- sprintf("%.17g",psi)
    if (exists(key,cache,inherits=FALSE)) {
      counts$candidate.cache.hits <<- counts$candidate.cache.hits+1L
      return(get(key,cache,inherits=FALSE))
    }
    z <- fit.profile(psi,y,"observed candidate")
    assign(key,z,cache); z
  }
  anchor <- candidate(psi.hat)
  evaluate.acceptability <- function(psi, B, phase) {
    check.deadline("acceptability")
    z <- candidate(psi)
    ans <- list(psi=psi,acceptability=NA_real_,status="observed_fit_failed",B=B,
                failed.draws=0L,statistic=NA_real_,phase=phase)
    if (isTRUE(z$converged) && isTRUE(anchor$converged)) {
      ans$statistic <- 2*(z$value-anchor$value)
      if (abs(psi-psi.hat) <= 1e-12*max(1,abs(psi.hat))) {
        ans$acceptability <- 1; ans$status <- "available"; ans$B <- 0L
      } else {
        rr <- .rd_exact_risks(drop(va%*%z$alpha),drop(vb%*%z$beta))
        prob <- pmin(pmax(ifelse(x==1,rr$p1,rr$p0),1e-12),1-1e-12)
        statistics <- rep(NA_real_,B)
        for (b in seq_len(B)) {
          check.deadline("bootstrap")
          yy <- rbinom(n,1,prob)
          a <- fit.profile(psi.hat,yy,"bootstrap fixed anchor")
          c <- fit.profile(psi,yy,"bootstrap candidate")
          counts$bootstrap.draws <<- counts$bootstrap.draws+1L
          if (isTRUE(a$converged) && isTRUE(c$converged)) statistics[b] <- 2*(c$value-a$value)
          else counts$failed.bootstrap.draws <<- counts$failed.bootstrap.draws+1L
        }
        ans$failed.draws <- sum(!is.finite(statistics))
        ## Never silently discard failed fits: that changes the null distribution.
        if (ans$failed.draws) ans$status <- "bootstrap_fit_failed" else {
          p.obs <- min(mean(statistics <= ans$statistic),mean(statistics >= ans$statistic))
          left <- rank(statistics,ties.method="max")/B
          right <- rank(-statistics,ties.method="max")/B
          ans$acceptability <- mean(pmin(left,right) <= p.obs)
          ans$status <- "available"
        }
      }
    }
    history[[length(history)+1L]] <<- as.data.frame(ans,stringsAsFactors=FALSE)
    ans
  }
  endpoint <- function(boundary) {
    if (!isTRUE(anchor$converged)) return(list(value=NA_real_,status="observed_fit_failed",available=FALSE))
    fractions <- .04
    while(tail(fractions,1) < 1) fractions <- c(fractions,min(1,tail(fractions,1)*grid.mult))
    grid <- psi.hat + fractions*(boundary-psi.hat)
    inner <- psi.hat; outer <- NA_real_; checked <- vector("list",length(grid))
    fail <- function(z) list(value=NA_real_,status=z$status,available=FALSE,
                            bracket=c(inner=inner,outer=outer),failure=z)
    for (k in seq_along(grid)) {
      z <- checked[[k]]
      if (is.null(z)) z <- evaluate.acceptability(grid[k],nsim.coarse,"coarse")
      if (!is.finite(z$acceptability)) return(fail(z))
      if (z$acceptability >= (1-level)-1e-12) { inner <- grid[k]; next }
      z <- evaluate.acceptability(grid[k],nsim.confirm,"confirm")
      if (!is.finite(z$acceptability)) return(fail(z))
      if (z$acceptability >= (1-level)-1e-12) { inner <- grid[k]; next }
      if (k < length(grid)) {
        next.z <- evaluate.acceptability(grid[k+1L],nsim.confirm,"confirm next")
        checked[[k+1L]] <- next.z
        if (!is.finite(next.z$acceptability)) return(fail(next.z))
        if (next.z$acceptability >= (1-level)-1e-12) next
      }
      outer <- grid[k]; break
    }
    if (!is.finite(outer)) return(list(value=boundary,status="boundary_truncated",available=TRUE,
      boundary.truncated=TRUE,converged=FALSE,bracket=c(inner=inner,outer=NA_real_),steps=0L))
    step <- 0L
    while(abs(inner-outer) > thres.dicho && step < bisection.max.step) {
      step <- step+1L; mid <- (inner+outer)/2
      z <- evaluate.acceptability(mid,nsim.refine,"refine")
      if (!is.finite(z$acceptability)) return(fail(z))
      if (z$acceptability >= (1-level)-1e-12) inner <- mid else outer <- mid
    }
    done <- abs(inner-outer) <= thres.dicho
    list(value=if(done) inner else NA_real_,last.accepted=inner,
         status=if(done) "available" else "search_not_converged",available=done,
         boundary.truncated=FALSE,converged=done,bracket=c(inner=inner,outer=outer),steps=step)
  }
  not.run <- list(value=NA_real_,status="not_requested",available=FALSE)
  lower <- upper <- not.run
  if (compute.interval) { lower <- endpoint(psi.range[1]); upper <- endpoint(psi.range[2]) }
  grid.results <- lapply(psi.grid,function(p) evaluate.acceptability(p,nsim.refine,"requested grid"))
  pv <- if(compute.p.value) evaluate.acceptability(psi.null,nsim.p.value,"p value") else
    list(acceptability=NA_real_,status="not_requested")
  interval.available <- isTRUE(lower$available) && isTRUE(upper$available)
  status <- if (!compute.interval) "not_requested" else if (!interval.available) "interval_failed" else
    if (isTRUE(lower$boundary.truncated) || isTRUE(upper$boundary.truncated)) "boundary_truncated" else "available"
  list(estimate=psi.hat,low=lower$value,up=upper$value,p=pv$acceptability,
       effect.estimate=psi.hat,effect.low=lower$value,effect.up=upper$value,
       level=level,target="marginal RD",interval="Blaker-type parametric-bootstrap inversion",
       conditioning="fixed covariates, treatment, fitting weights and standardization weights",
       status=status,interval.available=interval.available,p.status=pv$status,
       diagnostics=list(settings=list(optim.maxit=optim.maxit,optim.reltol=optim.reltol,
         optim.retries=optim.retries,optim.retry.maxit=optim.retry.maxit,
         optim.multistart=optim.multistart,optim.multistart.max=optim.multistart.max,
         optim.fallback=optim.fallback,optimizer.version="rd_bc_multistart_nlminb_v1",
         final.version="RD_final_B600_F1000_20261024",
         cutoff.rule="acceptability >= (1-level)-1e-12",cutoff.tolerance=1e-12,
         optim.boundary.rescue=optim.boundary.rescue,
         boundary.rescue.version="rd_bc_stable_boundary_v1",
         score.tol=score.tol,root.tol=root.tol,constraint.tol=constraint.tol,
         alpha.bound=alpha.bound,beta.bound=beta.bound,intercept=intercept,
         psi.range=psi.range,psi.se=psi.se,psi.se.source=psi.se.source,
         level=level,psi.null=psi.null,thres.dicho=thres.dicho,bisection.max.step=bisection.max.step,
         nsim=c(coarse=nsim.coarse,confirm=nsim.confirm,refine=nsim.refine,p.value=nsim.p.value),
         observed.anchor="fixed observed marginal estimate, including an FC input estimate",
         nuisance.fit="ordinary weighted likelihood; not a penalized FC likelihood",
         weight.uncertainty=FALSE,target.distribution.uncertainty=FALSE,
         search="outward grid and local inversion around anchor; disconnected sets not exhaustively searched"),
         target.at.input=target,anchor=anchor,lower=lower,upper=upper,p.value=pv,
         requested.grid=grid.results,candidates=as.list(cache),counts=counts,
         failed.fits=failures,retried.fits=retries,
         evaluations=if(length(history)) do.call(rbind,history) else data.frame(),
         elapsed.sec=as.numeric(Sys.time())-started))
}
