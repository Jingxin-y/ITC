## -------------------------------------------------------------------------
## Internal worker used by the optional bootstrap parallelization.
##
## This worker keeps the same bootstrap statistic and optimizer. Both the
## bootstrap null and alternative fits start from the original estimate
## passed into exact(); no warm start is used. It returns aggregate fitting
## diagnostics and does NOT retry failed fits (the original safe fallback is
## retained).
## -------------------------------------------------------------------------
.exact_bootstrap_lrt_worker <- function(y.local, state) {
  tryCatch(
    {
      check_deadline <- function(stage = "bootstrap_worker") {
        now <- as.numeric(Sys.time())
        if (now < state$exact.deadline) {
          return(invisible(NULL))
        }
        stop(
          sprintf(
            "exact() exceeded its %.0f-second time limit at stage '%s' (elapsed %.1f seconds)",
            state$time.limit.sec,
            stage,
            now - state$exact.start.time
          ),
          call. = FALSE
        )
      }

      diag_zero_worker <- function() {
        list(
          nuisance.fits = 0L,
          profile.steps = 0L,
          optim.calls = 0L,
          optim.function.evals = 0L,
          optimizer.failures = 0L,
          failed.nuisance.fits = 0L,
          retried.fits = 0L
        )
      }

      diag_add_worker <- function(x, y) {
        for (nm in names(x)) {
          x[[nm]] <- x[[nm]] + y[[nm]]
        }
        x
      }

      nll_fun_worker <- function(alpha, beta, yy) {
        p0p1 <- state$getProb(
          state$mat_vec_mul(state$va, alpha),
          state$mat_vec_mul(state$vb, beta)
        )

        p0 <- p0p1[state$idx0, 1]
        p1 <- p0p1[state$idx1, 2]

        p0 <- pmin(pmax(p0, state$eps), 1 - state$eps)
        p1 <- pmin(pmax(p1, state$eps), 1 - state$eps)

        y0.local <- yy[state$idx0]
        y1.local <- yy[state$idx1]

        -sum((1 - y0.local) * log(1 - p0) * state$w0 +
          y0.local * log(p0) * state$w0) -
          sum((1 - y1.local) * log(1 - p1) * state$w1 +
            y1.local * log(p1) * state$w1)
      }

      Diff_worker <- function(x1, x0) {
        sum((x1 - x0)^2) / sum(x1^2 + state$thres)
      }

      safe_optim_worker <- function(par, fn, bound) {
        check_deadline("bootstrap_safe_optim_before")

        lower <- rep(-bound, length(par))
        upper <- rep(bound, length(par))

        fit.raw <- tryCatch(
          stats::optim(
            par,
            fn,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(
              maxit = state$optim.maxit,
              factr = state$optim.reltol / .Machine$double.eps
            )
          ),
          error = function(e) {
            if (grepl("time limit", conditionMessage(e), ignore.case = TRUE)) {
              stop(e)
            }
            NULL
          }
        )

        check_deadline("bootstrap_safe_optim_after")

        fit.was.null <- is.null(fit.raw)
        nonfinite.par <- !fit.was.null && any(!is.finite(fit.raw$par))
        fallback.eval <- 0L

        ## Retain the original fallback behavior.  This is not a retry.
        if (fit.was.null || nonfinite.par) {
          fallback.value <- fn(par)
          fallback.eval <- 1L
          fit <- list(
            par = par,
            value = fallback.value,
            convergence = 99,
            counts = structure(c(0L, 0L), names = c("function", "gradient"))
          )
        } else {
          fit <- fit.raw
        }

        fn.evals <- 0L
        if (!is.null(fit.raw) && !is.null(fit.raw$counts)) {
          if (!is.null(names(fit.raw$counts)) && "function" %in% names(fit.raw$counts)) {
            fn.evals <- as.integer(fit.raw$counts[["function"]])
          } else if (length(fit.raw$counts) >= 1L) {
            fn.evals <- as.integer(fit.raw$counts[[1L]])
          }
        }
        if (is.na(fn.evals)) fn.evals <- 0L
        fn.evals <- fn.evals + fallback.eval

        convergence.bad <- is.null(fit$convergence) || !identical(as.integer(fit$convergence), 0L)
        value.bad <- is.null(fit$value) || length(fit$value) != 1L || !is.finite(fit$value)

        list(
          fit = fit,
          diagnostics = list(
            nuisance.fits = 0L,
            profile.steps = 0L,
            optim.calls = 1L,
            optim.function.evals = fn.evals,
            optimizer.failures = as.integer(
              fit.was.null || nonfinite.par || convergence.bad || value.bad
            ),
            failed.nuisance.fits = 0L,
            retried.fits = 0L
          )
        )
      }

      profile_fit_worker <- function(alphaj, j, yy, alpha.init, beta.init) {
        alpha <- alpha.init
        beta <- beta.init
        alpha[j] <- alphaj

        diff <- state$thres + 1
        step <- 0L
        diag <- diag_zero_worker()

        neg.log.likelihood.alpha <- function(alpha.in) {
          alpha.in[j] <- alphaj
          nll_fun_worker(alpha.in, beta, yy)
        }

        neg.log.likelihood.beta <- function(beta.in) {
          nll_fun_worker(alpha, beta.in, yy)
        }

        while (diff > state$thres && step < state$max.step) {
          check_deadline("bootstrap_profile_optimization")
          step <- step + 1L

          opt1.res <- safe_optim_worker(alpha, neg.log.likelihood.alpha, 8)
          diag <- diag_add_worker(diag, opt1.res$diagnostics)
          opt1 <- opt1.res$fit

          diff1 <- Diff_worker(opt1$par, alpha)
          alpha <- opt1$par
          alpha[j] <- alphaj

          opt2.res <- safe_optim_worker(beta, neg.log.likelihood.beta, 10)
          diag <- diag_add_worker(diag, opt2.res$diagnostics)
          opt2 <- opt2.res$fit

          diff2 <- Diff_worker(opt2$par, beta)
          beta <- opt2$par

          diff <- max(diff1, diff2)
        }

        converged <- is.finite(diff) && diff <= state$thres
        diag$nuisance.fits <- diag$nuisance.fits + 1L
        diag$profile.steps <- diag$profile.steps + step
        diag$failed.nuisance.fits <- diag$failed.nuisance.fits + as.integer(!converged)

        list(
          value = nll_fun_worker(alpha, beta, yy),
          alpha = alpha,
          beta = beta,
          converged = converged,
          profile.steps = step,
          diagnostics = diag
        )
      }

      ## Bootstrap null fit: start from the observed-data null constrained fit.
      null.fit <- profile_fit_worker(
        state$alpha.ml[state$j],
        state$j,
        y.local,
        alpha.init = state$alpha.null.start,
        beta.init = state$beta.null.start
      )

      ## Bootstrap alternative fit: use the same original-input start as the
      ## bootstrap null fit.
      alt.fit <- profile_fit_worker(
        state$alphaj,
        state$j,
        y.local,
        alpha.init = state$alpha.null.start,
        beta.init = state$beta.null.start
      )

      diag <- diag_add_worker(null.fit$diagnostics, alt.fit$diagnostics)
      value <- 2 * (alt.fit$value - null.fit$value)

      list(
        ok = TRUE,
        value = value,
        diagnostics = diag,
        message = NA_character_
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        value = NA_real_,
        diagnostics = NULL,
        message = conditionMessage(e)
      )
    }
  )
}

exact <- function(param, y, x, va, vb, weight = NULL,
                  max.step, thres = 1e-3, thres.dicho = 1e-3,
                  pars, se, pa, pb, optim.maxit = 300,
                  optim.reltol = 1e-6,
                  time.limit.sec = 10800,
                  parallel.bootstrap = FALSE,
                  ncores.bootstrap = 1L,
                  bisection.max.step = 15L,
                  grid.mult = 2) {
  time.limit.sec <- as.numeric(time.limit.sec)[1L]
  if (is.na(time.limit.sec) || time.limit.sec <= 0) {
    stop("time.limit.sec must be a positive number of seconds.", call. = FALSE)
  }

  bisection.max.step <- as.integer(bisection.max.step)[1L]
  if (is.na(bisection.max.step) || bisection.max.step < 1L) {
    stop("bisection.max.step must be a positive integer.", call. = FALSE)
  }

  grid.mult <- as.numeric(grid.mult)[1L]
  if (is.na(grid.mult) || !is.finite(grid.mult) || grid.mult <= 1) {
    stop("grid.mult must be a finite number greater than 1.", call. = FALSE)
  }

  thres.dicho <- as.numeric(thres.dicho)[1L]
  if (is.na(thres.dicho) || !is.finite(thres.dicho) || thres.dicho <= 0) {
    stop("thres.dicho must be a positive finite number.", call. = FALSE)
  }

  parallel.bootstrap <- isTRUE(parallel.bootstrap)
  ncores.bootstrap <- as.integer(ncores.bootstrap)[1L]
  if (is.na(ncores.bootstrap) || ncores.bootstrap < 1L) {
    stop("ncores.bootstrap must be a positive integer", call. = FALSE)
  }
  if (!parallel.bootstrap || ncores.bootstrap == 1L) {
    parallel.bootstrap <- FALSE
    ncores.bootstrap <- 1L
  }

  ## Cross-platform backend:
  ##   Windows       -> PSOCK cluster
  ##   Linux/macOS   -> fork (mclapply)
  ##   parallel=FALSE -> serial
  parallel.backend <- if (!parallel.bootstrap) {
    "serial"
  } else if (.Platform$OS.type == "windows") {
    "psock"
  } else {
    "fork"
  }

  ## Use wall-clock time so the same deadline can be checked inside
  ## independent PSOCK workers on Windows as well as forked workers on Unix.
  exact.start.time <- as.numeric(Sys.time())
  exact.deadline <- if (is.infinite(time.limit.sec)) {
    Inf
  } else {
    exact.start.time + time.limit.sec
  }

  check_exact_deadline <- function(stage = "unknown") {
    now <- as.numeric(Sys.time())
    if (now < exact.deadline) {
      return(invisible(NULL))
    }

    timeout.condition <- structure(
      list(
        message = sprintf(
          "exact() exceeded its %.0f-second time limit at stage '%s' (elapsed %.1f seconds)",
          time.limit.sec, stage, now - exact.start.time
        ),
        call = NULL,
        time.limit.sec = time.limit.sec,
        elapsed = now - exact.start.time,
        stage = stage
      ),
      class = c("exact_timeout", "error", "condition")
    )
    stop(timeout.condition)
  }

  if (is.null(weight)) {
    weight <- rep(1, length(y))
  }

  ## ------------------------------------------------------------
  ## Setup
  ## ------------------------------------------------------------
  getProb <- if (param == "RR") getProbRR else getProbRD

  ## Create the Windows PSOCK cluster only once per exact() call and reuse it
  ## for every coarse/confirm/refine bootstrap evaluation.
  bootstrap.cl <- NULL
  if (parallel.backend == "psock") {
    bootstrap.cl <- parallel::makeCluster(ncores.bootstrap)
    on.exit(
      {
        try(parallel::stopCluster(bootstrap.cl), silent = TRUE)
      },
      add = TRUE
    )

    ## If mat_vec_mul() is provided by a package/Rcpp DLL, load that package
    ## on each clean PSOCK worker before exporting the wrapper.  The package
    ## name is detected from the function environment or from a PACKAGE= entry
    ## in the generated .Call wrapper.
    matmul.pkg <- NULL
    if (exists("mat_vec_mul", mode = "function", inherits = TRUE)) {
      matmul.fun.master <- get("mat_vec_mul", mode = "function", inherits = TRUE)
      env.name <- environmentName(environment(matmul.fun.master))
      if (startsWith(env.name, "namespace:")) {
        matmul.pkg <- sub("^namespace:", "", env.name)
      } else {
        body.txt <- paste(deparse(body(matmul.fun.master)), collapse = " ")
        mm <- regexec(
          "PACKAGE\\s*=\\s*['\"]([^'\"]+)['\"]",
          body.txt,
          perl = TRUE
        )
        hit <- regmatches(body.txt, mm)[[1L]]
        if (length(hit) >= 2L) matmul.pkg <- hit[2L]
      }
    }

    if (is.null(matmul.pkg) && "brmplus" %in% loadedNamespaces()) {
      matmul.pkg <- "brmplus"
    }

    if (!is.null(matmul.pkg) && nzchar(matmul.pkg)) {
      pkg.ok <- parallel::clusterCall(
        bootstrap.cl,
        function(pkg) {
          suppressPackageStartupMessages(
            require(pkg, character.only = TRUE, quietly = TRUE)
          )
        },
        matmul.pkg
      )
      if (!all(vapply(pkg.ok, isTRUE, logical(1)))) {
        stop(
          "Could not load package '", matmul.pkg,
          "' on all Windows PSOCK workers.",
          call. = FALSE
        )
      }
    }

    ## PSOCK workers start from clean global environments.  Export the helper
    ## functions used by the sourced brmplus code whenever they exist in the
    ## master's global environment.  Missing names are simply ignored.
    helper.names <- c(
      "getProbRR", "getProbRD",
      "getProbScalarRR", "getProbScalarRD",
      "same", "getPrbAux", "mat_vec_mul"
    )
    helper.names <- helper.names[
      vapply(
        helper.names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
      )
    ]
    if (length(helper.names)) {
      parallel::clusterExport(
        bootstrap.cl,
        varlist = helper.names,
        envir = .GlobalEnv
      )
    }
  }

  ## Build a self-contained copy of getProb for PSOCK workers.  Functions
  ## sourced into .GlobalEnv (notably same() and getProbScalar*) are copied
  ## into a small private environment, so the worker does not depend on the
  ## master's global environment layout.
  make_portable_getProb <- function(fun) {
    helper.env <- new.env(parent = baseenv())
    helper.names.local <- c(
      "same", "getPrbAux", "getProbScalarRR", "getProbScalarRD"
    )

    source.env <- environment(fun)
    for (nm in helper.names.local) {
      if (exists(nm, envir = source.env, inherits = TRUE)) {
        obj <- get(nm, envir = source.env, inherits = TRUE)
        if (is.function(obj) && identical(typeof(obj), "closure")) {
          obj <- unserialize(serialize(obj, NULL))
          environment(obj) <- helper.env
        }
        assign(nm, obj, envir = helper.env)
      }
    }

    out <- unserialize(serialize(fun, NULL))
    if (is.function(out) && identical(typeof(out), "closure")) {
      environment(out) <- helper.env
    }
    out
  }

  getProb.worker <- if (parallel.backend == "psock") {
    make_portable_getProb(getProb)
  } else {
    getProb
  }

  ## Keep the same matrix-vector multiplication routine used by the original
  ## exact() implementation.  The design matrices themselves are created once
  ## and reused; only parameter-dependent products are recomputed.
  mat_vec_mul.worker <- mat_vec_mul

  alpha.ml <- pars[1:pa]
  beta.ml <- pars[(pa + 1):(pa + pb)]

  ## Precompute indices and invariant subsets once.
  idx0 <- x == 0
  idx1 <- x == 1

  y0 <- y[idx0]
  y1 <- y[idx1]

  w0 <- weight[idx0]
  w1 <- weight[idx1]

  va0 <- va[idx0, , drop = FALSE]
  va1 <- va[idx1, , drop = FALSE]
  vb0 <- vb[idx0, , drop = FALSE]
  vb1 <- vb[idx1, , drop = FALSE]

  n0 <- sum(idx0)
  n1 <- sum(idx1)
  n <- length(y)

  eps <- 1e-12

  ## Baseline start used only to obtain the observed-data null constrained fit.
  alpha.start <- alpha.ml
  beta.start <- beta.ml

  ## ------------------------------------------------------------
  ## Diagnostics helpers
  ## ------------------------------------------------------------
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
      retried.fits = 0L,
      candidate.fits = 0L,
      candidate.cache.hits = 0L,
      acceptability.evaluations = 0L,
      grid.evaluations = 0L,
      bisection.evaluations = 0L,
      bisection.steps = 0L
    )
  }

  diag_add <- function(x, y) {
    for (nm in names(x)) {
      if (!is.null(y[[nm]])) {
        x[[nm]] <- x[[nm]] + y[[nm]]
      }
    }
    x
  }

  ## Convert worker diagnostics to the full master diagnostic shape.
  expand_worker_diag <- function(d) {
    out <- diag_zero()
    if (is.null(d)) {
      return(out)
    }
    for (nm in intersect(names(d), names(out))) {
      out[[nm]] <- d[[nm]]
    }
    out
  }

  ## ------------------------------------------------------------
  ## Negative log-likelihood helper
  ## ------------------------------------------------------------
  nll_fun <- function(alpha, beta, y.local) {
    p0p1 <- getProb(
      mat_vec_mul(va, alpha),
      mat_vec_mul(vb, beta)
    )

    p0 <- p0p1[idx0, 1]
    p1 <- p0p1[idx1, 2]

    p0 <- pmin(pmax(p0, eps), 1 - eps)
    p1 <- pmin(pmax(p1, eps), 1 - eps)

    y0.local <- y.local[idx0]
    y1.local <- y.local[idx1]

    -sum((1 - y0.local) * log(1 - p0) * w0 +
      y0.local * log(p0) * w0) -
      sum((1 - y1.local) * log(1 - p1) * w1 +
        y1.local * log(p1) * w1)
  }

  Diff <- function(x1, x0) {
    sum((x1 - x0)^2) / sum(x1^2 + thres)
  }

  safe_optim <- function(par, fn, bound) {
    check_exact_deadline("safe_optim_before")
    lower <- rep(-bound, length(par))
    upper <- rep(bound, length(par))

    fit.raw <- tryCatch(
      stats::optim(
        par,
        fn,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = optim.maxit,
          factr = optim.reltol / .Machine$double.eps
        )
      ),
      error = function(e) {
        ## Do not swallow elapsed-time errors.  Ordinary optimization failures
        ## retain the original fallback below; there is deliberately no retry.
        if (inherits(e, "exact_timeout") ||
          grepl("time limit", conditionMessage(e), ignore.case = TRUE)) {
          stop(e)
        }
        NULL
      }
    )
    check_exact_deadline("safe_optim_after")

    fit.was.null <- is.null(fit.raw)
    nonfinite.par <- !fit.was.null && any(!is.finite(fit.raw$par))
    fallback.eval <- 0L

    ## Original fallback behavior.  This is NOT a failure retry.
    if (fit.was.null || nonfinite.par) {
      fallback.value <- fn(par)
      fallback.eval <- 1L
      fit <- list(
        par = par,
        value = fallback.value,
        convergence = 99,
        counts = structure(c(0L, 0L), names = c("function", "gradient"))
      )
    } else {
      fit <- fit.raw
    }

    fn.evals <- 0L
    if (!is.null(fit.raw) && !is.null(fit.raw$counts)) {
      if (!is.null(names(fit.raw$counts)) && "function" %in% names(fit.raw$counts)) {
        fn.evals <- as.integer(fit.raw$counts[["function"]])
      } else if (length(fit.raw$counts) >= 1L) {
        fn.evals <- as.integer(fit.raw$counts[[1L]])
      }
    }
    if (is.na(fn.evals)) fn.evals <- 0L
    fn.evals <- fn.evals + fallback.eval

    convergence.bad <- is.null(fit$convergence) || !identical(as.integer(fit$convergence), 0L)
    value.bad <- is.null(fit$value) || length(fit$value) != 1L || !is.finite(fit$value)

    out.diag <- diag_zero()
    out.diag$optim.calls <- 1L
    out.diag$optim.function.evals <- fn.evals
    out.diag$optimizer.failures <- as.integer(
      fit.was.null || nonfinite.par || convergence.bad || value.bad
    )

    list(fit = fit, diagnostics = out.diag)
  }

  ## ------------------------------------------------------------
  ## Unified profile nuisance optimization
  ## ------------------------------------------------------------
  ## This replaces the duplicated observed-data constrained fit that used to
  ## appear once in LRT.alpha() and again at the start of ptail().
  profile_fit <- function(alphaj, j, y.local,
                          alpha.init = alpha.start,
                          beta.init = beta.start,
                          stage = "profile_optimization") {
    alpha <- alpha.init
    beta <- beta.init
    alpha[j] <- alphaj

    diff <- thres + 1
    step <- 0L
    fit.diag <- diag_zero()

    neg.log.likelihood.alpha <- function(alpha.in) {
      alpha.in[j] <- alphaj
      nll_fun(alpha.in, beta, y.local)
    }

    neg.log.likelihood.beta <- function(beta.in) {
      nll_fun(alpha, beta.in, y.local)
    }

    while (diff > thres && step < max.step) {
      check_exact_deadline(stage)
      step <- step + 1L

      opt1.res <- safe_optim(alpha, neg.log.likelihood.alpha, 8)
      fit.diag <- diag_add(fit.diag, opt1.res$diagnostics)
      opt1 <- opt1.res$fit

      diff1 <- Diff(opt1$par, alpha)
      alpha <- opt1$par
      alpha[j] <- alphaj

      opt2.res <- safe_optim(beta, neg.log.likelihood.beta, 10)
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
      alpha = alpha,
      beta = beta,
      converged = converged,
      profile.steps = step,
      diagnostics = fit.diag
    )
  }

  ## ------------------------------------------------------------
  ## Cache the observed-data null profile once per alpha component
  ## ------------------------------------------------------------
  ## This fit is the fixed center anchor for all observed-data candidate fits
  ## and the matched starting value for bootstrap null fits.
  null.fit.obs <- vector("list", pa)
  null.fit.diag <- vector("list", pa)

  for (j in seq_len(pa)) {
    check_exact_deadline("observed_null_profile")
    null.fit.obs[[j]] <- profile_fit(
      alpha.ml[j],
      j,
      y,
      alpha.init = alpha.start,
      beta.init = beta.start,
      stage = "observed_null_profile"
    )
    null.fit.diag[[j]] <- null.fit.obs[[j]]$diagnostics
  }

  ## ------------------------------------------------------------
  ## Simulate distribution of profile-LRT statistic
  ## ------------------------------------------------------------
  ## candidate.fit is the already-computed observed constrained fit at alphaj.
  ## Hence ptail() no longer refits the same observed nuisance parameters.
  ptail <- function(alphaj, j, candidate.fit, nsim = 500) {
    ## Fitted probabilities under the observed constrained candidate fit.
    prob <- getProb(
      mat_vec_mul(va, candidate.fit$alpha),
      mat_vec_mul(vb, candidate.fit$beta)
    )

    p0 <- prob[idx0, 1]
    p1 <- prob[idx1, 2]

    p0 <- pmin(pmax(p0, eps), 1 - eps)
    p1 <- pmin(pmax(p1, eps), 1 - eps)

    ## Generate all bootstrap outcomes in the parent process first.  This keeps
    ## serial and parallel runs on the same RNG stream for a fixed call path.
    y.sim.list <- vector("list", nsim)

    for (i in seq_len(nsim)) {
      if (i %% 10L == 1L) {
        check_exact_deadline("parametric_bootstrap_generation")
      }

      y.sim <- numeric(n)
      y.sim[idx0] <- rbinom(n0, 1, p0)
      y.sim[idx1] <- rbinom(n1, 1, p1)
      y.sim.list[[i]] <- y.sim
    }

    boot.diag <- diag_zero()
    boot.diag$bootstrap.draws <- as.integer(nsim)

    if (parallel.backend == "serial") {
      one_bootstrap_lrt <- function(y.sim) {
        ## Null fit starts from the original estimate passed into exact().
        null.sim <- profile_fit(
          alpha.ml[j],
          j,
          y.sim,
          alpha.init = alpha.start,
          beta.init = beta.start,
          stage = "bootstrap_null_profile"
        )

        ## Alternative fit uses the same original-input start as the null fit.
        alt.sim <- profile_fit(
          alphaj,
          j,
          y.sim,
          alpha.init = alpha.start,
          beta.init = beta.start,
          stage = "bootstrap_alt_profile"
        )

        list(
          value = 2 * (alt.sim$value - null.sim$value),
          diagnostics = diag_add(null.sim$diagnostics, alt.sim$diagnostics)
        )
      }

      ans <- lapply(y.sim.list, one_bootstrap_lrt)
      LRT.sim <- vapply(ans, `[[`, numeric(1), "value")
      for (ii in seq_along(ans)) {
        boot.diag <- diag_add(boot.diag, ans[[ii]]$diagnostics)
      }
    } else {
      worker.state <- list(
        getProb = getProb.worker,
        mat_vec_mul = mat_vec_mul.worker,
        va = va,
        vb = vb,
        idx0 = idx0,
        idx1 = idx1,
        w0 = w0,
        w1 = w1,
        eps = eps,
        thres = thres,
        max.step = max.step,
        optim.maxit = optim.maxit,
        optim.reltol = optim.reltol,
        alpha.null.start = alpha.start,
        beta.null.start = beta.start,
        alpha.ml = alpha.ml,
        alphaj = alphaj,
        j = j,
        exact.start.time = exact.start.time,
        exact.deadline = exact.deadline,
        time.limit.sec = time.limit.sec
      )

      if (parallel.backend == "fork") {
        ans <- parallel::mclapply(
          y.sim.list,
          .exact_bootstrap_lrt_worker,
          state = worker.state,
          mc.cores = min(ncores.bootstrap, nsim),
          mc.preschedule = TRUE,
          mc.set.seed = FALSE
        )
      } else {
        ans <- parallel::parLapply(
          bootstrap.cl,
          y.sim.list,
          .exact_bootstrap_lrt_worker,
          state = worker.state
        )
      }

      ok.boot <- vapply(ans, `[[`, logical(1), "ok")
      if (any(!ok.boot)) {
        msg <- paste(
          unique(vapply(ans[!ok.boot], `[[`, character(1), "message")),
          collapse = " | "
        )

        if (grepl("time limit", msg, ignore.case = TRUE)) {
          now <- as.numeric(Sys.time())
          timeout.condition <- structure(
            list(
              message = msg,
              call = NULL,
              time.limit.sec = time.limit.sec,
              elapsed = now - exact.start.time,
              stage = paste0("parametric_bootstrap_", parallel.backend)
            ),
            class = c("exact_timeout", "error", "condition")
          )
          stop(timeout.condition)
        }

        stop(
          "bootstrap LRT calculation failed on ", parallel.backend,
          " backend: ", msg,
          call. = FALSE
        )
      }

      LRT.sim <- vapply(ans, `[[`, numeric(1), "value")
      for (ii in seq_along(ans)) {
        boot.diag <- diag_add(
          boot.diag,
          expand_worker_diag(ans[[ii]]$diagnostics)
        )
      }
    }

    check_exact_deadline("parametric_bootstrap_complete")
    list(values = LRT.sim, diagnostics = boot.diag)
  }

  ## ------------------------------------------------------------
  ## Acceptability function (unchanged statistical definition)
  ## ------------------------------------------------------------
  acceptability <- function(alphaj, LRT.obs, LRT.sim) {
    check_exact_deadline("acceptability")

    p.left.obs <- mean(LRT.sim <= LRT.obs)
    p.right.obs <- mean(LRT.sim >= LRT.obs)
    p.min.obs <- min(p.left.obs, p.right.obs)

    p.left <- vapply(LRT.sim, function(z) mean(LRT.sim <= z), numeric(1))
    p.right <- vapply(LRT.sim, function(z) mean(LRT.sim >= z), numeric(1))
    p.min <- pmin(p.left, p.right)

    mean(p.min <= p.min.obs)
  }

  ## ------------------------------------------------------------
  ## Dichotomy with 2x outward grid, fixed original-input candidate starts,
  ## cached observed candidate fits, and explicit bisection diagnostics.
  ## ------------------------------------------------------------
  dichotomy <- function(j, alpha.low, alpha.up,
                        direction = "low",
                        thres.dicho.local = thres.dicho,
                        max.bisection = bisection.max.step,
                        grid.mult.local = grid.mult,
                        nsim.coarse = 100,
                        nsim.confirm = 200,
                        nsim.refine = 200) {
    endpoint.diag <- diag_zero()

    mle <- if (direction == "low") alpha.up else alpha.low
    boundary <- if (direction == "low") alpha.low else alpha.up

    ## Candidate-fit cache. Seed it with the observed-data null fit at the
    ## MLE-side point. New candidate fits do NOT inherit another candidate's
    ## nuisance solution; every new candidate starts from this same center fit.
    candidate.values <- mle
    candidate.fits <- list(null.fit.obs[[j]])

    cache_tol <- function(a) {
      1e-12 * max(1, abs(a))
    }

    find_exact_candidate <- function(a) {
      if (!length(candidate.values)) {
        return(NA_integer_)
      }
      d <- abs(candidate.values - a)
      k <- which.min(d)
      if (length(k) && d[k] <= cache_tol(a)) k else NA_integer_
    }

    get_candidate_fit <- function(a) {
      cached <- find_exact_candidate(a)
      if (!is.na(cached)) {
        endpoint.diag$candidate.cache.hits <<-
          endpoint.diag$candidate.cache.hits + 1L
        return(candidate.fits[[cached]])
      }

      ## Fixed original-input start: every observed-data candidate is fit
      ## independently from the estimate originally passed into exact().
      fit <- profile_fit(
        a,
        j,
        y,
        alpha.init = alpha.start,
        beta.init = beta.start,
        stage = "observed_candidate_profile"
      )

      endpoint.diag <<- diag_add(endpoint.diag, fit$diagnostics)
      endpoint.diag$candidate.fits <<-
        endpoint.diag$candidate.fits + 1L

      candidate.values <<- c(candidate.values, a)
      candidate.fits[[length(candidate.fits) + 1L]] <<- fit
      fit
    }

    eval_accept <- function(a, nsim, phase = c("coarse", "confirm", "refine")) {
      phase <- match.arg(phase)
      check_exact_deadline("grid_or_refinement_evaluation")

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
      LRT.obs <- 2 * (fit.a$value - null.fit.obs[[j]]$value)

      boot <- ptail(
        a,
        j,
        candidate.fit = fit.a,
        nsim = nsim
      )
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

      acceptability(a, LRT.obs, boot$values)
    }

    ## Search fractions are now exactly 0, .04, .08, .16, .32, .64, 1 when
    ## grid.mult.local = 2.
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

    ## MLE-side point.
    a.vals[1] <- eval_accept(grid[1], nsim = nsim.coarse, phase = "coarse")
    if (a.vals[1] < 0.05) {
      a.vals[1] <- eval_accept(grid[1], nsim = nsim.confirm, phase = "confirm")
    }

    if (a.vals[1] < 0.05) {
      endpoint.diag$bisection.steps <- 0L
      return(list(
        alpha.dicho = grid[1],
        convergence = FALSE,
        reason = "MLE-side point rejected",
        grid = grid,
        a.vals = a.vals,
        bracket = NULL,
        candidate.values = candidate.values,
        candidate.converged = vapply(candidate.fits, function(z) isTRUE(z$converged), logical(1)),
        diagnostics = endpoint.diag
      ))
    }

    inner.idx <- 1L
    outer.idx <- NA_integer_

    for (k in 2:length(grid)) {
      ## A point may already have been evaluated as the outward confirmation
      ## point in the preceding iteration.
      if (is.na(a.vals[k])) {
        a.vals[k] <- eval_accept(grid[k], nsim = nsim.coarse, phase = "coarse")
      }

      if (a.vals[k] >= 0.05) {
        inner.idx <- k
        next
      }

      ## Re-evaluate a coarse rejection with B=200 confirmation.
      a.vals[k] <- eval_accept(grid[k], nsim = nsim.confirm, phase = "confirm")
      if (a.vals[k] >= 0.05) {
        inner.idx <- k
        next
      }

      ## Keep the existing two-consecutive-rejection stability rule.
      if (k < length(grid)) {
        a.next <- eval_accept(
          grid[k + 1L],
          nsim = nsim.confirm,
          phase = "confirm"
        )
        a.vals[k + 1L] <- a.next
        if (a.next < 0.05) {
          outer.idx <- k
          break
        }
        inner.idx <- k + 1L
      } else {
        outer.idx <- k
        break
      }
    }

    ## If no confirmed rejection was found, the endpoint is outside the range.
    if (is.na(outer.idx)) {
      evaluated <- which(!is.na(a.vals))
      endpoint.diag$bisection.steps <- 0L
      return(list(
        alpha.dicho = grid[length(grid)],
        convergence = FALSE,
        reason = "confirmed rejected point not found in search range",
        grid = grid[evaluated],
        a.vals = a.vals[evaluated],
        bracket = NULL,
        candidate.values = candidate.values,
        candidate.converged = vapply(candidate.fits, function(z) isTRUE(z$converged), logical(1)),
        diagnostics = endpoint.diag
      ))
    }

    ## Bracket: inner accepted, outer rejected.
    inner <- grid[inner.idx]
    outer <- grid[outer.idx]

    step <- 0L

    ## Stop once width <= thres.dicho.local or after max.bisection steps.
    while (abs(inner - outer) > thres.dicho.local &&
      step < max.bisection) {
      check_exact_deadline("dichotomy")

      step <- step + 1L
      mid <- (inner + outer) / 2

      a.mid <- eval_accept(mid, nsim = nsim.refine, phase = "refine")

      if (a.mid >= 0.05) {
        inner <- mid
      } else {
        outer <- mid
      }
    }

    endpoint.diag$bisection.steps <- step

    evaluated <- which(!is.na(a.vals))
    list(
      alpha.dicho = inner,
      convergence = abs(inner - outer) <= thres.dicho.local,
      reason = if (abs(inner - outer) <= thres.dicho.local) {
        NA_character_
      } else {
        "maximum bisection steps reached before width threshold"
      },
      grid = grid[evaluated],
      a.vals = a.vals[evaluated],
      bracket = c(inner = inner, outer = outer),
      candidate.values = candidate.values,
      candidate.converged = vapply(candidate.fits, function(z) isTRUE(z$converged), logical(1)),
      diagnostics = endpoint.diag
    )
  }

  ## ------------------------------------------------------------
  ## Get candidate alpha bounds
  ## ------------------------------------------------------------
  cilen <- ifelse(length(y) < 300, 8, 4)
  alpha.up.start <- pmin(alpha.ml + cilen * se[1:pa], 8)
  alpha.low.start <- pmax(alpha.ml - cilen * se[1:pa], -8)

  alpha.up1 <- rep(0, pa)
  alpha.low1 <- rep(0, pa)
  up.res <- vector("list", pa)
  low.res <- vector("list", pa)

  for (j in seq_len(pa)) {
    check_exact_deadline("upper_endpoint")
    up.res[[j]] <- dichotomy(
      j,
      alpha.ml[j],
      alpha.up.start[j],
      direction = "up",
      thres.dicho.local = thres.dicho,
      max.bisection = bisection.max.step,
      grid.mult.local = grid.mult
    )
    alpha.up1[j] <- up.res[[j]]$alpha.dicho
  }

  for (j in seq_len(pa)) {
    check_exact_deadline("lower_endpoint")
    low.res[[j]] <- dichotomy(
      j,
      alpha.low.start[j],
      alpha.ml[j],
      direction = "low",
      thres.dicho.local = thres.dicho,
      max.bisection = bisection.max.step,
      grid.mult.local = grid.mult
    )
    alpha.low1[j] <- low.res[[j]]$alpha.dicho
  }

  ## ------------------------------------------------------------
  ## P-values: same B=500 and same acceptability definition.
  ## ------------------------------------------------------------
  p.value <- rep(0, pa)
  p.diag <- vector("list", pa)

  for (j in seq_len(pa)) {
    check_exact_deadline("p_value")
    this.diag <- diag_zero()

    ## Fit the observed alpha=0 profile independently from the original
    ## estimate passed into exact(). This keeps the observed P-value fit on
    ## the same fixed-start rule as all non-bootstrap candidate fits.
    fit.p <- profile_fit(
      0,
      j,
      y,
      alpha.init = alpha.start,
      beta.init = beta.start,
      stage = "p_value_observed_profile"
    )
    this.diag <- diag_add(this.diag, fit.p$diagnostics)
    this.diag$candidate.fits <- this.diag$candidate.fits + 1L

    LRT.obs.p <- 2 * (fit.p$value - null.fit.obs[[j]]$value)

    boot.p <- ptail(
      0,
      j,
      candidate.fit = fit.p,
      nsim = 500
    )
    this.diag <- diag_add(this.diag, boot.p$diagnostics)
    this.diag$bootstrap.draws.p.value <-
      this.diag$bootstrap.draws.p.value + 500L
    this.diag$acceptability.evaluations <-
      this.diag$acceptability.evaluations + 1L

    p.value[j] <- acceptability(
      0,
      LRT.obs.p,
      boot.p$values
    )

    p.diag[[j]] <- this.diag
  }

  ## ------------------------------------------------------------
  ## Aggregate diagnostics without changing the old low/up/p interface.
  ## ------------------------------------------------------------
  total.diag <- diag_zero()

  for (j in seq_len(pa)) {
    total.diag <- diag_add(total.diag, null.fit.diag[[j]])
    total.diag <- diag_add(total.diag, up.res[[j]]$diagnostics)
    total.diag <- diag_add(total.diag, low.res[[j]]$diagnostics)
    total.diag <- diag_add(total.diag, p.diag[[j]])
  }

  diagnostics <- list(
    settings = list(
      thres = thres,
      thres.dicho = thres.dicho,
      bisection.max.step = bisection.max.step,
      grid.mult = grid.mult,
      grid.fractions = {
        ff <- c(0, 0.04)
        while (tail(ff, 1) < 1) {
          ff <- c(ff, min(1, tail(ff, 1) * grid.mult))
        }
        unique(ff)
      },
      nsim.coarse = 100L,
      nsim.confirm = 200L,
      nsim.refine = 200L,
      nsim.p.value = 500L,
      parallel.backend = parallel.backend,
      ncores.bootstrap = ncores.bootstrap,
      failure.retry.enabled = FALSE,
      observed.candidate.start = "original_input_estimate",
      bootstrap.null.start = "original_input_estimate",
      bootstrap.candidate.start = "original_input_estimate",
      pvalue.observed.start = "original_input_estimate",
      pvalue.bootstrap.null.start = "original_input_estimate",
      pvalue.bootstrap.candidate.start = "original_input_estimate"
    ),
    observed.null = lapply(
      seq_len(pa),
      function(j) {
        list(
          j = j,
          value = null.fit.obs[[j]]$value,
          converged = null.fit.obs[[j]]$converged,
          diagnostics = null.fit.diag[[j]]
        )
      }
    ),
    upper = lapply(
      seq_len(pa),
      function(j) {
        list(
          j = j,
          endpoint = up.res[[j]]$alpha.dicho,
          convergence = up.res[[j]]$convergence,
          reason = up.res[[j]]$reason,
          candidate.converged = up.res[[j]]$candidate.converged,
          bracket = up.res[[j]]$bracket,
          final.bracket.width = if (is.null(up.res[[j]]$bracket)) {
            NA_real_
          } else {
            abs(diff(up.res[[j]]$bracket))
          },
          diagnostics = up.res[[j]]$diagnostics
        )
      }
    ),
    lower = lapply(
      seq_len(pa),
      function(j) {
        list(
          j = j,
          endpoint = low.res[[j]]$alpha.dicho,
          convergence = low.res[[j]]$convergence,
          reason = low.res[[j]]$reason,
          candidate.converged = low.res[[j]]$candidate.converged,
          bracket = low.res[[j]]$bracket,
          final.bracket.width = if (is.null(low.res[[j]]$bracket)) {
            NA_real_
          } else {
            abs(diff(low.res[[j]]$bracket))
          },
          diagnostics = low.res[[j]]$diagnostics
        )
      }
    ),
    p.value = lapply(
      seq_len(pa),
      function(j) {
        list(
          j = j,
          p = p.value[j],
          diagnostics = p.diag[[j]]
        )
      }
    ),
    total = total.diag,
    elapsed.sec = as.numeric(Sys.time()) - exact.start.time
  )

  list(
    low = alpha.low1,
    up = alpha.up1,
    p = p.value,
    diagnostics = diagnostics
  )
}
