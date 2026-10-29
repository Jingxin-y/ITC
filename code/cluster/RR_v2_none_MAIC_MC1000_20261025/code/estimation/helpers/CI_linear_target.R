## -------------------------------------------------------------------------
## Blaker-type parametric-bootstrap CI for a linear target
##
##   psi = contrast^T alpha
##
## This generalizes CI_exact_adapt_para.R, where the constrained target is
## one coordinate alpha_j.  Here the constrained target can be any nonzero
## linear combination of alpha.
##
## The statistical construction is intentionally kept parallel to exact():
##   1) For each candidate psi0, profile all nuisance parameters under
##          contrast^T alpha = psi0.
##   2) Generate bootstrap outcomes from that candidate-specific constrained
##      fit.
##   3) In each bootstrap sample, refit both the fixed observed-data anchor
##      psi_hat = contrast^T alpha_hat and the candidate psi0.
##   4) Use the same Blaker-type acceptability calculation.
##   5) Invert acceptability >= level using outward search + bisection.
##
## Required brmplus helpers in the calling environment:
##   getProbRR(), getProbRD()
## Optional:
##   mat_vec_mul(); if unavailable, base R %*% is used.
## -------------------------------------------------------------------------

exact_linear_target <- function(
    param,
    y,
    x,
    va,
    vb,
    contrast,
    weight = NULL,
    max.step,
    pars,
    pa,
    pb,
    vcov.alpha = NULL,
    psi.se = NULL,
    se = NULL,
    level = 0.95,
    psi.null = 0,
    thres = 1e-3,
    thres.dicho = 1e-3,
    optim.maxit = 50,
    optim.reltol = 1e-6,
    time.limit.sec = 10800,
    bisection.max.step = 15L,
    grid.mult = 2,
    nsim.coarse = 100L,
    nsim.confirm = 200L,
    nsim.refine = 200L,
    nsim.p.value = 500L,
    alpha.bound = 8,
    beta.bound = 10,
    gamma.bound = NULL,
    search.mult = NULL,
    seed = NULL
) {
  ## -----------------------------------------------------------------------
  ## Input checks
  ## -----------------------------------------------------------------------
  param <- match.arg(param, c("RR", "RD"))

  y <- as.numeric(y)
  x <- as.numeric(x)
  va <- as.matrix(va)
  vb <- as.matrix(vb)
  contrast <- as.numeric(contrast)
  pars <- as.numeric(pars)

  if (length(contrast) != pa) {
    stop("length(contrast) must equal pa.", call. = FALSE)
  }
  if (any(!is.finite(contrast)) || sum(contrast^2) <= 0) {
    stop("contrast must be a finite, nonzero vector.", call. = FALSE)
  }
  if (length(pars) < pa + pb) {
    stop("pars must contain pa alpha parameters followed by pb beta parameters.",
         call. = FALSE)
  }
  if (nrow(va) != length(y) || nrow(vb) != length(y) || length(x) != length(y)) {
    stop("y, x, va, and vb must have compatible dimensions.", call. = FALSE)
  }
  if (ncol(va) != pa || ncol(vb) != pb) {
    stop("ncol(va) must equal pa and ncol(vb) must equal pb.", call. = FALSE)
  }

  if (is.null(weight)) weight <- rep(1, length(y))
  weight <- as.numeric(weight)
  if (length(weight) != length(y)) {
    stop("weight must have the same length as y.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  alpha.level <- 1 - level
  if (!is.finite(alpha.level) || alpha.level <= 0 || alpha.level >= 1) {
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  }

  if (is.null(search.mult)) search.mult <- ifelse(length(y) < 300, 8, 4)
  if (is.null(gamma.bound)) gamma.bound <- alpha.bound * sqrt(max(pa, 1))

  ## -----------------------------------------------------------------------
  ## Deadline
  ## -----------------------------------------------------------------------
  exact.start.time <- as.numeric(Sys.time())
  exact.deadline <- if (is.infinite(time.limit.sec)) {
    Inf
  } else {
    exact.start.time + as.numeric(time.limit.sec)
  }

  check_deadline <- function(stage = "unknown") {
    now <- as.numeric(Sys.time())
    if (now < exact.deadline) return(invisible(NULL))
    stop(
      sprintf(
        "exact_linear_target() exceeded its %.0f-second time limit at stage '%s' (elapsed %.1f seconds)",
        time.limit.sec, stage, now - exact.start.time
      ),
      call. = FALSE
    )
  }

  ## -----------------------------------------------------------------------
  ## Model setup
  ## -----------------------------------------------------------------------
  getProb <- if (param == "RR") getProbRR else getProbRD

  mmul <- if (exists("mat_vec_mul", mode = "function", inherits = TRUE)) {
    get("mat_vec_mul", mode = "function", inherits = TRUE)
  } else {
    function(A, b) drop(A %*% b)
  }

  alpha.ml <- pars[seq_len(pa)]
  beta.ml <- pars[pa + seq_len(pb)]

  psi.hat <- sum(contrast * alpha.ml)

  idx0 <- x == 0
  idx1 <- x == 1
  n0 <- sum(idx0)
  n1 <- sum(idx1)
  n <- length(y)
  w0 <- weight[idx0]
  w1 <- weight[idx1]
  eps <- 1e-12

  alpha.start <- alpha.ml
  beta.start <- beta.ml

  ## -----------------------------------------------------------------------
  ## Reparameterization under contrast^T alpha = psi
  ##
  ## alpha = alpha_particular(psi) + N gamma,
  ## where columns of N form an orthonormal basis for null(contrast^T).
  ## This uses all alpha coefficients jointly and therefore preserves their
  ## dependence through the likelihood rather than treating them separately.
  ## -----------------------------------------------------------------------
  contrast.norm2 <- sum(contrast^2)

  alpha_particular <- function(psi) {
    contrast * (psi / contrast.norm2)
  }

  if (pa == 1L) {
    N <- matrix(numeric(0), nrow = 1L, ncol = 0L)
  } else {
    ## First column of Q spans contrast; remaining columns span its orthogonal
    ## complement.  cbind(..., diag(pa)) guarantees a full p-dimensional QR.
    Q <- qr.Q(qr(cbind(contrast, diag(pa))), complete = TRUE)
    N <- Q[, -1L, drop = FALSE]
  }

  alpha_from_gamma <- function(psi, gamma) {
    a0 <- alpha_particular(psi)
    if (ncol(N) == 0L) return(as.numeric(a0))
    as.numeric(a0 + N %*% gamma)
  }

  gamma_from_alpha <- function(alpha, psi) {
    if (ncol(N) == 0L) return(numeric(0))
    as.numeric(crossprod(N, alpha - alpha_particular(psi)))
  }

  ## -----------------------------------------------------------------------
  ## Diagnostics
  ## -----------------------------------------------------------------------
  diag_zero <- function() {
    list(
      bootstrap.draws = 0L,
      bootstrap.draws.coarse = 0L,
      bootstrap.draws.confirm = 0L,
      bootstrap.draws.refine = 0L,
      bootstrap.draws.p.value = 0L,
      nuisance.fits = 0L,
      profile.steps = 0L,
      optim.calls = 0L,
      optim.function.evals = 0L,
      optimizer.failures = 0L,
      failed.nuisance.fits = 0L,
      candidate.fits = 0L,
      candidate.cache.hits = 0L,
      acceptability.evaluations = 0L,
      grid.evaluations = 0L,
      bisection.evaluations = 0L,
      bisection.steps = 0L
    )
  }

  diag_add <- function(x, z) {
    for (nm in intersect(names(x), names(z))) x[[nm]] <- x[[nm]] + z[[nm]]
    x
  }

  ## -----------------------------------------------------------------------
  ## Negative log likelihood
  ## -----------------------------------------------------------------------
  nll_fun <- function(alpha, beta, y.local) {
    p0p1 <- getProb(mmul(va, alpha), mmul(vb, beta))

    p0 <- pmin(pmax(p0p1[idx0, 1], eps), 1 - eps)
    p1 <- pmin(pmax(p0p1[idx1, 2], eps), 1 - eps)

    y0.local <- y.local[idx0]
    y1.local <- y.local[idx1]

    -sum((1 - y0.local) * log(1 - p0) * w0 + y0.local * log(p0) * w0) -
      sum((1 - y1.local) * log(1 - p1) * w1 + y1.local * log(p1) * w1)
  }

  Diff <- function(x1, x0) {
    if (!length(x1)) return(0)
    sum((x1 - x0)^2) / sum(x1^2 + thres)
  }

  safe_optim <- function(par, fn, bound) {
    check_deadline("safe_optim_before")

    if (!length(par)) {
      value <- fn(par)
      d <- diag_zero()
      return(list(
        fit = list(par = par, value = value, convergence = 0L,
                   counts = c("function" = 1L, gradient = 0L)),
        diagnostics = d
      ))
    }

    fit.raw <- tryCatch(
      stats::optim(
        par,
        fn,
        method = "L-BFGS-B",
        lower = rep(-bound, length(par)),
        upper = rep(bound, length(par)),
        control = list(
          maxit = optim.maxit,
          factr = optim.reltol / .Machine$double.eps
        )
      ),
      error = function(e) NULL
    )

    d <- diag_zero()
    d$optim.calls <- 1L

    if (is.null(fit.raw) || any(!is.finite(fit.raw$par))) {
      fallback.value <- fn(par)
      fit <- list(
        par = par,
        value = fallback.value,
        convergence = 99L,
        counts = c("function" = 1L, gradient = 0L)
      )
      d$optim.function.evals <- 1L
      d$optimizer.failures <- 1L
    } else {
      fit <- fit.raw
      fn.evals <- if (!is.null(fit$counts)) as.integer(fit$counts[[1L]]) else 0L
      if (is.na(fn.evals)) fn.evals <- 0L
      d$optim.function.evals <- fn.evals
      if (!identical(as.integer(fit$convergence), 0L) || !is.finite(fit$value)) {
        d$optimizer.failures <- 1L
      }
    }

    check_deadline("safe_optim_after")
    list(fit = fit, diagnostics = d)
  }

  ## -----------------------------------------------------------------------
  ## Constrained profile likelihood for a candidate psi0
  ##
  ## Free nuisance parameters:
  ##   gamma : pa - 1 directions orthogonal to contrast
  ##   beta  : GOP nuisance coefficients
  ## -----------------------------------------------------------------------
  profile_fit_target <- function(
      psi0,
      y.local,
      alpha.init = alpha.start,
      beta.init = beta.start,
      stage = "profile_target"
  ) {
    gamma <- gamma_from_alpha(alpha.init, psi0)
    beta <- beta.init
    alpha <- alpha_from_gamma(psi0, gamma)

    diff <- thres + 1
    step <- 0L
    fit.diag <- diag_zero()

    while (diff > thres && step < max.step) {
      check_deadline(stage)
      step <- step + 1L

      ## Optimize alpha nuisance directions while enforcing c^T alpha = psi0.
      if (length(gamma)) {
        neg.log.likelihood.gamma <- function(gamma.in) {
          alpha.in <- alpha_from_gamma(psi0, gamma.in)
          nll_fun(alpha.in, beta, y.local)
        }

        opt1.res <- safe_optim(gamma, neg.log.likelihood.gamma, gamma.bound)
        fit.diag <- diag_add(fit.diag, opt1.res$diagnostics)
        opt1 <- opt1.res$fit

        diff1 <- Diff(opt1$par, gamma)
        gamma <- opt1$par
        alpha <- alpha_from_gamma(psi0, gamma)
      } else {
        diff1 <- 0
        alpha <- alpha_from_gamma(psi0, gamma)
      }

      ## Optimize GOP nuisance parameters.
      neg.log.likelihood.beta <- function(beta.in) {
        nll_fun(alpha, beta.in, y.local)
      }

      opt2.res <- safe_optim(beta, neg.log.likelihood.beta, beta.bound)
      fit.diag <- diag_add(fit.diag, opt2.res$diagnostics)
      opt2 <- opt2.res$fit

      diff2 <- Diff(opt2$par, beta)
      beta <- opt2$par

      diff <- max(diff1, diff2)
    }

    converged <- is.finite(diff) && diff <= thres
    fit.diag$nuisance.fits <- fit.diag$nuisance.fits + 1L
    fit.diag$profile.steps <- fit.diag$profile.steps + step
    fit.diag$failed.nuisance.fits <-
      fit.diag$failed.nuisance.fits + as.integer(!converged)

    list(
      value = nll_fun(alpha, beta, y.local),
      psi = psi0,
      alpha = alpha,
      gamma = gamma,
      beta = beta,
      constraint.error = sum(contrast * alpha) - psi0,
      converged = converged,
      profile.steps = step,
      diagnostics = fit.diag
    )
  }

  ## -----------------------------------------------------------------------
  ## Fixed observed-data anchor
  ## -----------------------------------------------------------------------
  check_deadline("observed_anchor_profile")
  anchor.fit.obs <- profile_fit_target(
    psi.hat,
    y,
    alpha.init = alpha.start,
    beta.init = beta.start,
    stage = "observed_anchor_profile"
  )

  ## -----------------------------------------------------------------------
  ## Candidate-specific parametric bootstrap distribution
  ## -----------------------------------------------------------------------
  ptail_target <- function(psi0, candidate.fit, nsim = 500L) {
    prob <- getProb(
      mmul(va, candidate.fit$alpha),
      mmul(vb, candidate.fit$beta)
    )

    p0 <- pmin(pmax(prob[idx0, 1], eps), 1 - eps)
    p1 <- pmin(pmax(prob[idx1, 2], eps), 1 - eps)

    LRT.sim <- rep(NA_real_, nsim)
    boot.diag <- diag_zero()
    boot.diag$bootstrap.draws <- as.integer(nsim)

    for (i in seq_len(nsim)) {
      if (i %% 10L == 1L) check_deadline("parametric_bootstrap")

      y.sim <- numeric(n)
      y.sim[idx0] <- stats::rbinom(n0, 1, p0)
      y.sim[idx1] <- stats::rbinom(n1, 1, p1)

      ## Fixed observed-data anchor, exactly analogous to the current code.
      anchor.sim <- profile_fit_target(
        psi.hat,
        y.sim,
        alpha.init = alpha.start,
        beta.init = beta.start,
        stage = "bootstrap_anchor_profile"
      )

      candidate.sim <- profile_fit_target(
        psi0,
        y.sim,
        alpha.init = alpha.start,
        beta.init = beta.start,
        stage = "bootstrap_candidate_profile"
      )

      boot.diag <- diag_add(boot.diag, anchor.sim$diagnostics)
      boot.diag <- diag_add(boot.diag, candidate.sim$diagnostics)

      LRT.sim[i] <- 2 * (candidate.sim$value - anchor.sim$value)
    }

    list(values = LRT.sim, diagnostics = boot.diag)
  }

  ## -----------------------------------------------------------------------
  ## Blaker-type acceptability -- same definition as current exact()
  ## -----------------------------------------------------------------------
  acceptability <- function(LRT.obs, LRT.sim) {
    check_deadline("acceptability")

    p.left.obs <- mean(LRT.sim <= LRT.obs)
    p.right.obs <- mean(LRT.sim >= LRT.obs)
    p.min.obs <- min(p.left.obs, p.right.obs)

    p.left <- vapply(LRT.sim, function(z) mean(LRT.sim <= z), numeric(1))
    p.right <- vapply(LRT.sim, function(z) mean(LRT.sim >= z), numeric(1))
    p.min <- pmin(p.left, p.right)

    mean(p.min <= p.min.obs)
  }

  ## -----------------------------------------------------------------------
  ## One-dimensional inversion in psi
  ## -----------------------------------------------------------------------
  dichotomy_target <- function(
      psi.low,
      psi.up,
      direction = c("low", "up"),
      thres.dicho.local = thres.dicho,
      max.bisection = bisection.max.step,
      grid.mult.local = grid.mult
  ) {
    direction <- match.arg(direction)
    endpoint.diag <- diag_zero()

    mle <- if (direction == "low") psi.up else psi.low
    boundary <- if (direction == "low") psi.low else psi.up

    candidate.values <- mle
    candidate.fits <- list(anchor.fit.obs)

    cache_tol <- function(a) 1e-12 * max(1, abs(a))

    find_candidate <- function(a) {
      d <- abs(candidate.values - a)
      k <- which.min(d)
      if (length(k) && d[k] <= cache_tol(a)) k else NA_integer_
    }

    get_candidate_fit <- function(a) {
      cached <- find_candidate(a)
      if (!is.na(cached)) {
        endpoint.diag$candidate.cache.hits <<-
          endpoint.diag$candidate.cache.hits + 1L
        return(candidate.fits[[cached]])
      }

      fit <- profile_fit_target(
        a,
        y,
        alpha.init = alpha.start,
        beta.init = beta.start,
        stage = "observed_candidate_profile"
      )

      endpoint.diag <<- diag_add(endpoint.diag, fit$diagnostics)
      endpoint.diag$candidate.fits <<- endpoint.diag$candidate.fits + 1L
      candidate.values <<- c(candidate.values, a)
      candidate.fits[[length(candidate.fits) + 1L]] <<- fit
      fit
    }

    eval_accept <- function(a, nsim, phase = c("coarse", "confirm", "refine")) {
      phase <- match.arg(phase)
      check_deadline("grid_or_refinement_evaluation")

      endpoint.diag$acceptability.evaluations <<-
        endpoint.diag$acceptability.evaluations + 1L
      if (phase == "refine") {
        endpoint.diag$bisection.evaluations <<-
          endpoint.diag$bisection.evaluations + 1L
      } else {
        endpoint.diag$grid.evaluations <<-
          endpoint.diag$grid.evaluations + 1L
      }

      fit.a <- get_candidate_fit(a)
      LRT.obs <- 2 * (fit.a$value - anchor.fit.obs$value)

      boot <- ptail_target(a, candidate.fit = fit.a, nsim = nsim)
      endpoint.diag <<- diag_add(endpoint.diag, boot$diagnostics)

      if (phase == "coarse") {
        endpoint.diag$bootstrap.draws.coarse <<-
          endpoint.diag$bootstrap.draws.coarse + as.integer(nsim)
      } else if (phase == "confirm") {
        endpoint.diag$bootstrap.draws.confirm <<-
          endpoint.diag$bootstrap.draws.confirm + as.integer(nsim)
      } else {
        endpoint.diag$bootstrap.draws.refine <<-
          endpoint.diag$bootstrap.draws.refine + as.integer(nsim)
      }

      acceptability(LRT.obs, boot$values)
    }

    fractions <- c(0, 0.04)
    while (tail(fractions, 1) < 1) {
      fractions <- c(
        fractions,
        min(1, tail(fractions, 1) * grid.mult.local)
      )
    }
    fractions <- unique(fractions)

    grid <- mle + fractions * (boundary - mle)
    a.vals <- rep(NA_real_, length(grid))

    ## Anchor-side point should be accepted.
    a.vals[1] <- eval_accept(grid[1], nsim.coarse, "coarse")
    if (a.vals[1] < alpha.level) {
      a.vals[1] <- eval_accept(grid[1], nsim.confirm, "confirm")
    }

    if (a.vals[1] < alpha.level) {
      return(list(
        psi.dicho = grid[1],
        convergence = FALSE,
        reason = "anchor-side point rejected",
        grid = grid[1],
        a.vals = a.vals[1],
        bracket = NULL,
        diagnostics = endpoint.diag
      ))
    }

    inner.idx <- 1L
    outer.idx <- NA_integer_

    for (k in 2:length(grid)) {
      if (is.na(a.vals[k])) {
        a.vals[k] <- eval_accept(grid[k], nsim.coarse, "coarse")
      }

      if (a.vals[k] >= alpha.level) {
        inner.idx <- k
        next
      }

      ## Confirm coarse rejection.
      a.vals[k] <- eval_accept(grid[k], nsim.confirm, "confirm")
      if (a.vals[k] >= alpha.level) {
        inner.idx <- k
        next
      }

      ## Preserve the current two-consecutive-rejection rule.
      if (k < length(grid)) {
        a.next <- eval_accept(grid[k + 1L], nsim.confirm, "confirm")
        a.vals[k + 1L] <- a.next
        if (a.next < alpha.level) {
          outer.idx <- k
          break
        }
        inner.idx <- k + 1L
      } else {
        outer.idx <- k
        break
      }
    }

    if (is.na(outer.idx)) {
      evaluated <- which(!is.na(a.vals))
      return(list(
        psi.dicho = grid[length(grid)],
        convergence = FALSE,
        reason = "confirmed rejected point not found in search range",
        grid = grid[evaluated],
        a.vals = a.vals[evaluated],
        bracket = NULL,
        diagnostics = endpoint.diag
      ))
    }

    inner <- grid[inner.idx]
    outer <- grid[outer.idx]
    step <- 0L

    while (abs(inner - outer) > thres.dicho.local && step < max.bisection) {
      check_deadline("dichotomy")
      step <- step + 1L
      mid <- (inner + outer) / 2
      a.mid <- eval_accept(mid, nsim.refine, "refine")

      if (a.mid >= alpha.level) {
        inner <- mid
      } else {
        outer <- mid
      }
    }

    endpoint.diag$bisection.steps <- step
    evaluated <- which(!is.na(a.vals))

    list(
      psi.dicho = inner,
      convergence = abs(inner - outer) <= thres.dicho.local,
      reason = if (abs(inner - outer) <= thres.dicho.local) NA_character_ else
        "maximum bisection steps reached before width threshold",
      grid = grid[evaluated],
      a.vals = a.vals[evaluated],
      bracket = c(inner = inner, outer = outer),
      diagnostics = endpoint.diag
    )
  }

  ## -----------------------------------------------------------------------
  ## Starting search width
  ##
  ## This SE is used ONLY to define the numerical outward-search range.
  ## The reported CI itself is obtained by bootstrap inversion.
  ## -----------------------------------------------------------------------
  psi.se.source <- NULL

  if (is.null(psi.se) && !is.null(vcov.alpha)) {
    vcov.alpha <- as.matrix(vcov.alpha)
    if (!all(dim(vcov.alpha) == c(pa, pa))) {
      stop("vcov.alpha must be a pa x pa covariance matrix.", call. = FALSE)
    }
    psi.var <- as.numeric(crossprod(contrast, vcov.alpha %*% contrast))
    psi.se <- sqrt(max(psi.var, 0))
    psi.se.source <- "contrast^T vcov.alpha contrast"
  } else if (!is.null(psi.se)) {
    psi.se <- as.numeric(psi.se)[1L]
    psi.se.source <- "user supplied"
  } else if (!is.null(se)) {
    ## Fallback only.  Correlations are ignored here, but this quantity affects
    ## only the initial search range, not the bootstrap-inversion definition.
    se <- as.numeric(se)
    if (length(se) < pa) stop("se must contain at least pa entries.", call. = FALSE)
    psi.se <- sqrt(sum((contrast * se[seq_len(pa)])^2))
    psi.se.source <- "diagonal SE fallback (search range only)"
    warning(
      "vcov.alpha was not supplied. psi.se was approximated from marginal SEs; ",
      "this affects only the initial search range. Supplying vcov.alpha is preferred.",
      call. = FALSE
    )
  } else {
    stop("Supply vcov.alpha, psi.se, or se to set the initial search range.",
         call. = FALSE)
  }

  if (!is.finite(psi.se) || psi.se <= 0) {
    stop("psi.se must be positive and finite.", call. = FALSE)
  }

  ## If each alpha coefficient is numerically bounded near +/- alpha.bound,
  ## the largest absolute linear target within that box is this value.
  psi.abs.max <- alpha.bound * sum(abs(contrast))

  psi.low.start <- max(psi.hat - search.mult * psi.se, -psi.abs.max)
  psi.up.start <- min(psi.hat + search.mult * psi.se, psi.abs.max)

  ## -----------------------------------------------------------------------
  ## CI endpoints
  ## -----------------------------------------------------------------------
  low.res <- dichotomy_target(
    psi.low = psi.low.start,
    psi.up = psi.hat,
    direction = "low"
  )

  up.res <- dichotomy_target(
    psi.low = psi.hat,
    psi.up = psi.up.start,
    direction = "up"
  )

  ## -----------------------------------------------------------------------
  ## P-value for H0: psi = psi.null
  ## -----------------------------------------------------------------------
  p.diag <- diag_zero()

  fit.p <- profile_fit_target(
    psi.null,
    y,
    alpha.init = alpha.start,
    beta.init = beta.start,
    stage = "p_value_observed_profile"
  )
  p.diag <- diag_add(p.diag, fit.p$diagnostics)
  p.diag$candidate.fits <- p.diag$candidate.fits + 1L

  LRT.obs.p <- 2 * (fit.p$value - anchor.fit.obs$value)
  boot.p <- ptail_target(psi.null, fit.p, nsim.p.value)
  p.diag <- diag_add(p.diag, boot.p$diagnostics)
  p.diag$bootstrap.draws.p.value <-
    p.diag$bootstrap.draws.p.value + as.integer(nsim.p.value)
  p.diag$acceptability.evaluations <- p.diag$acceptability.evaluations + 1L

  p.value <- acceptability(LRT.obs.p, boot.p$values)

  ## -----------------------------------------------------------------------
  ## Transform link-scale target to RR or RD scale
  ## -----------------------------------------------------------------------
  transform_effect <- if (param == "RR") exp else tanh

  diagnostics <- list(
    settings = list(
      level = level,
      contrast = contrast,
      psi.null = psi.null,
      thres = thres,
      thres.dicho = thres.dicho,
      bisection.max.step = bisection.max.step,
      grid.mult = grid.mult,
      nsim.coarse = nsim.coarse,
      nsim.confirm = nsim.confirm,
      nsim.refine = nsim.refine,
      nsim.p.value = nsim.p.value,
      search.mult = search.mult,
      psi.se = psi.se,
      psi.se.source = psi.se.source,
      alpha.bound = alpha.bound,
      beta.bound = beta.bound,
      gamma.bound = gamma.bound,
      observed.anchor = "psi_hat = contrast^T alpha_hat",
      candidate.constraint = "contrast^T alpha = psi0"
    ),
    anchor = list(
      psi = psi.hat,
      alpha = anchor.fit.obs$alpha,
      beta = anchor.fit.obs$beta,
      value = anchor.fit.obs$value,
      constraint.error = anchor.fit.obs$constraint.error,
      converged = anchor.fit.obs$converged,
      diagnostics = anchor.fit.obs$diagnostics
    ),
    lower = low.res,
    upper = up.res,
    p.value = list(
      psi.null = psi.null,
      p = p.value,
      candidate.fit = fit.p,
      diagnostics = p.diag
    ),
    elapsed.sec = as.numeric(Sys.time()) - exact.start.time
  )

  list(
    estimate = psi.hat,
    low = low.res$psi.dicho,
    up = up.res$psi.dicho,
    p = p.value,
    effect.estimate = transform_effect(psi.hat),
    effect.low = transform_effect(low.res$psi.dicho),
    effect.up = transform_effect(up.res$psi.dicho),
    contrast = contrast,
    diagnostics = diagnostics
  )
}

## -------------------------------------------------------------------------
## Examples of defining the target
## -------------------------------------------------------------------------
## If alpha(v0) = alpha_1 + v0 * alpha_2:
##
##   v0 <- 0.30
##   contrast <- c(1, v0)
##
##   fit.ci <- exact_linear_target(
##     param = "RR",
##     y = y,
##     x = x,
##     va = va,
##     vb = vb,
##     contrast = contrast,
##     weight = weight,
##     max.step = 50,
##     pars = pars,
##     pa = 2,
##     pb = ncol(vb),
##     vcov.alpha = V.alpha,
##     seed = 1
##   )
##
## fit.ci[c("estimate", "low", "up", "p")]
## fit.ci[c("effect.estimate", "effect.low", "effect.up")]
##
## Sanity check against the old coordinate-wise implementation:
##   contrast = c(1, 0, ..., 0)
## should give approximately the old CI for alpha_1, up to Monte Carlo noise
## and small numerical differences from the more explicit reparameterization.
