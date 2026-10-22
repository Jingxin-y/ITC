## Shared RR/RD root solver for the fixed-analytic-weight adjusted score.
## Convergence requires a small score residual, not merely a small parameter step.
## Armijo uses the actual search direction. Oversized/ill-conditioned Newton
## steps use an SVD damped least-squares direction, without changing the score.
max.likelihood.firth <- function(y, x, va, vb, alpha.start, beta.start, weight,
                                 max.step, thres, pa, pb, derivatives) {
  n <- length(y)
  if (n == 0L || length(x) != n || length(weight) != n ||
      nrow(va) != n || nrow(vb) != n || ncol(va) != pa || ncol(vb) != pb ||
      any(!is.finite(c(y, x, va, vb, weight))) ||
      any(!y %in% c(0, 1)) || any(!x %in% c(0, 1)) ||
      any(weight < 0) || !any(weight > 0)) {
    stop("FC requires finite binary outcomes/treatments, conformable designs and nonnegative weights with positive sum.")
  }
  if (length(max.step) != 1L || !is.finite(max.step) || max.step < 1 ||
      length(thres) != 1L || !is.finite(thres) || thres <= 0) {
    stop("FC requires max.step >= 1 and a positive finite score tolerance.")
  }
  theta <- c(alpha.start, beta.start)
  if (length(alpha.start) != pa || length(beta.start) != pb ||
      any(!is.finite(theta))) stop("Invalid FC starting values.")

  ## Both U and A are homogeneous of degree one in the weights. Normalize
  ## internally so root tolerances and the numerical path share this invariance.
  weight <- weight / max(weight)
  weight <- weight / mean(weight)
  evaluations <- 0L
  evaluate <- function(par) {
    evaluations <<- evaluations + 1L
    terms <- derivatives(par[seq_len(pa)], par[pa + seq_len(pb)], x, va, vb)
    if (isFALSE(terms$interior) || any(!is.finite(terms$p)) ||
        any(terms$p <= 1e-8 | terms$p >= 1 - 1e-8)) {
      stop("FC candidate is outside the untruncated interior probability range.")
    }
    ans <- firth_score_adjustment(y, va, vb, weight, terms)
    if (any(!is.finite(c(ans$adjusted.score, ans$information, ans$covariance)))) {
      stop("Non-finite FC score or covariance.")
    }
    ans
  }
  invalid.candidates <- 0L
  last.error <- NULL
  candidate <- function(par) tryCatch(evaluate(par), error = function(e) {
    invalid.candidates <<- invalid.candidates + 1L
    last.error <<- conditionMessage(e)
    NULL
  })
  current <- evaluate(theta)
  score.scale <- pmax(1, sqrt(diag(current$information)))
  working.scale <- score.scale
  residual <- function(fit) fit$adjusted.score / working.scale
  root.residual <- function(fit) fit$adjusted.score / score.scale
  r <- residual(current)
  step <- 0L
  reason <- "maximum iterations reached"
  direction.counts <- c(newton=0L, damped=0L, gradient=0L)
  line.searches <- rejected.steps <- 0L
  max.condition <- 0
  last.stationarity <- NA_real_
  last.search <- NULL
  restarts <- 0L
  restart.reason <- NULL
  restart.starts <- NULL
  restart.index <- 0L
  restart.history <- list()
  attempt.start.step <- 0L
  attempt.limit <- Inf
  attempt.name <- "supplied"
  attempt.recorded <- FALSE
  best <- list(par=theta,fit=current,
               score=max(abs(root.residual(current))),working.scale=working.scale)
  remember <- function(par,fit) {
    score <- max(abs(root.residual(fit)))
    if (score < best$score) best <<- list(par=par,fit=fit,score=score,
                                        working.scale=working.scale)
  }
  record.attempt <- function(why) {
    restart.history[[length(restart.history)+1L]] <<- list(
      start=attempt.name,
      iterations=step-attempt.start.step,
      max_scaled_score=max(abs(root.residual(current))),reason=why)
    attempt.recorded <<- TRUE
  }

  ## Deterministic restarts after a stalled search. At alpha=0,
  ## both RR and RD mappings have p0=p1, hence log(OP)=2*logit(p).
  ## This changes only the initial values, not the FC score or fitting weights.
  ## The original score scale and TOTAL iteration budget remain fixed.
  ## The existing OP-only attempt comes first. Additional symmetric treatment-
  ## effect starts may cross a nonzero stationary point of the residual loss.
  ## Each additional start receives at most 50 of the remaining total iterations.
  make.starts <- function() {
    cell <- match(apply(vb,1,paste,collapse="/"),unique(apply(vb,1,paste,collapse="/")))
    probability <- vapply(seq_len(max(cell)),function(k) {
      ix <- cell==k & weight>0
      (sum(weight[ix]*y[ix])+.5)/(sum(weight[ix])+1)
    },numeric(1))
    logop <- 2*(log(probability[cell])-log1p(-probability[cell]))
    beta <- drop(qr.solve(vb*sqrt(weight),logop*sqrt(weight)))
    starts <- list(list(name="pooled_op",par=c(rep(0,pa),beta)))
    intercept <- which(vapply(seq_len(pa),function(j)
      all(abs(va[,j]-1)<1e-12),logical(1)))
    intercept <- if(length(intercept)) intercept[1L] else NA_integer_
    slopes <- setdiff(seq_len(pa),intercept[!is.na(intercept)])
    pairs <- matrix(c(-4,-4,4,4,0,-2,0,2,-2,-4,2,4),ncol=2,byrow=TRUE)
    if(length(slopes)) {
      for(j in slopes) for(k in seq_len(nrow(pairs))) {
        alpha <- rep(0,pa)
        if(!is.na(intercept)) alpha[intercept] <- pairs[k,1L]
        alpha[j] <- pairs[k,2L]
        starts[[length(starts)+1L]] <- list(
          name=paste0("effect_",j,"_",pairs[k,1L],"_",pairs[k,2L]),
          par=c(alpha,beta))
      }
    } else if(!is.na(intercept)) {
      for(value in c(-4,4,-2,2)) {
        alpha <- rep(0,pa); alpha[intercept] <- value
        starts[[length(starts)+1L]] <- list(
          name=paste0("intercept_",value),par=c(alpha,beta))
      }
    }
    starts
  }
  restart <- function(why) {
    if (step >= max.step) return(FALSE)
    record.attempt(why)
    if(is.null(restart.starts))
      restart.starts <<- tryCatch(make.starts(),error=function(e) list())
    while(restart.index < length(restart.starts)) {
      restart.index <<- restart.index+1L
      start <- restart.starts[[restart.index]]$par
      restarts <<- restarts+1L
      restart.reason <<- why
      fit <- if(all(is.finite(start))) candidate(start) else NULL
      if(is.null(fit)) {
        restart.history[[length(restart.history)+1L]] <<- list(
          start=restart.starts[[restart.index]]$name,iterations=0L,
          max_scaled_score=NA_real_,reason="invalid restart candidate")
        next
      }
      theta <<- start
      current <<- fit
      ## Rebalance only the merit function; final root acceptance always uses
      ## the ORIGINAL score scale, independently of the starting values.
      working.scale <<- pmax(1,sqrt(diag(fit$information)))
      r <<- residual(fit)
      remember(theta,current)
      attempt.start.step <<- step
      attempt.limit <<- if(restart.index==1L) Inf else 50L
      attempt.name <<- restart.starts[[restart.index]]$name
      attempt.recorded <<- FALSE
      return(TRUE)
    }
    FALSE
  }

  ## L = ||r||^2. The derivative of L along the ACTUAL (possibly capped)
  ## direction is 2*r'J*d; a fixed fraction of L is incorrect after capping.
  search <- function(direction, jacobian, loss, kind, damping=0) {
    if (any(!is.finite(direction)) || max(abs(direction)) == 0) return(NULL)
    direction <- direction / max(1, max(abs(direction))/2)
    slope <- 2 * sum(r * drop(jacobian %*% direction))
    if (!is.finite(slope) || slope >= 0) return(NULL)
    line.searches <<- line.searches + 1L
    for (halving in 0:24) {
      fraction <- 2^(-halving)
      next.theta <- theta + fraction * direction
      next.fit <- candidate(next.theta)
      if (is.null(next.fit)) next
      next.r <- residual(next.fit)
      next.loss <- sum(next.r^2)
      last.search <<- list(direction=kind, damping=damping, fraction=fraction,
        slope=slope, old.loss=loss, new.loss=next.loss,
        required.loss=loss+1e-4*fraction*slope,
        max.scaled.score=max(abs(root.residual(next.fit))))
      ## A rounded/nondecreasing objective must not rescue an unsolved root.
      if ((next.loss < loss && next.loss <= loss+1e-4*fraction*slope) ||
          max(abs(root.residual(next.fit))) <= thres) {
        return(list(par=next.theta, fit=next.fit, residual=next.r,
                    kind=kind, fraction=fraction, loss=next.loss, search=last.search))
      }
      rejected.steps <<- rejected.steps + 1L
    }
    NULL
  }

  while (max(abs(root.residual(current))) > thres && step < max.step) {
    if(step-attempt.start.step >= attempt.limit) {
      reason <- "additional start reached its iteration allowance"
      if(restart(reason)) { reason <- "maximum iterations reached"; next }
      break
    }
    step <- step + 1L
    jacobian <- matrix(NA_real_, length(theta), length(theta))
    for (j in seq_along(theta)) {
      h <- 1e-5 * max(1, abs(theta[j]))
      ## Retry a smaller difference only when both sides leave the domain.
      for (difference.retry in 0:4) {
        plus <- minus <- theta
        plus[j] <- plus[j] + h
        minus[j] <- minus[j] - h
        fp <- candidate(plus)
        fm <- candidate(minus)
        if (!is.null(fp) || !is.null(fm)) break
        h <- h/10
      }
      if (!is.null(fp) && !is.null(fm)) {
        jacobian[, j] <- (residual(fp) - residual(fm)) / (2 * h)
      } else if (!is.null(fp)) {
        jacobian[, j] <- (residual(fp) - r) / h
      } else if (!is.null(fm)) {
        jacobian[, j] <- (r - residual(fm)) / h
      }
    }
    if (any(!is.finite(jacobian))) {
      reason <- "could not differentiate the adjusted score in the interior"
      if (restart(reason)) { reason <- "maximum iterations reached"; next }
      break
    }
    decomposition <- svd(jacobian)
    largest <- max(decomposition$d)
    smallest <- min(decomposition$d)
    condition <- if (smallest > 0) largest/smallest else Inf
    max.condition <- max(max.condition, condition)
    loss <- sum(r^2)
    last.stationarity <- max(abs(drop(crossprod(jacobian,r)))) /
      max(largest*sqrt(loss),.Machine$double.xmin)
    ## A nearly singular Jacobian and an almost orthogonal residual indicate
    ## stagnation of the loss away from a root. Restart rather than spending
    ## the remaining budget on floating-point-sized decreases. This never
    ## declares convergence: the original score test remains mandatory.
    if(condition > 1e8 && is.finite(last.stationarity) &&
       last.stationarity <= sqrt(.Machine$double.eps)) {
      reason <- "adjusted-score residual is stationary away from zero"
      if(restart(reason)) { reason <- "maximum iterations reached"; next }
      break
    }
    accepted <- NULL
    ## A very large Newton direction is dominated by weak numerical directions.
    ## Prefer damping to shrinking every coordinate by the same huge factor.
    if (is.finite(condition) && condition <= 1e8) {
      direction <- tryCatch(as.vector(solve(jacobian, -r)), error=function(e) NULL)
      if (!is.null(direction) && all(is.finite(direction)) && max(abs(direction)) <= 2)
        accepted <- search(direction, jacobian, loss, "newton")
    }
    if ((is.null(accepted) || accepted$fraction < .25) && is.finite(largest) && largest > 0) {
      projected <- drop(crossprod(decomposition$u, r))
      damped <- function(mu) -drop(decomposition$v %*%
        (decomposition$d/(decomposition$d^2+mu) * projected))
      damping <- max(largest^2*1e-8, .Machine$double.xmin)
      direction <- damped(damping)
      for (increase in 0:24) {
        if (all(is.finite(direction)) && max(abs(direction)) <= 2) break
        damping <- damping*10
        direction <- damped(damping)
      }
      ## Try increasingly conservative directions on the SAME data and score.
      for (retry in 0:8) {
        proposal <- search(damped(damping), jacobian, loss, "damped", damping)
        if (!is.null(proposal) && (is.null(accepted) || proposal$loss < accepted$loss))
          accepted <- proposal
        ## A tiny accepted fraction can stagnate for hundreds of iterations.
        ## Try stronger damping; retain the best actual improvement found.
        if (!is.null(proposal) && (proposal$fraction >= .25 ||
                                  max(abs(root.residual(proposal$fit))) <= thres)) break
        damping <- damping*10
      }
    }
    if (is.null(accepted))
      accepted <- search(-as.vector(crossprod(jacobian,r)), jacobian, loss, "gradient")
    if (is.null(accepted)) {
      reason <- "line search did not reduce the adjusted-score residual"
      if (restart(reason)) { reason <- "maximum iterations reached"; next }
      break
    }
    theta <- accepted$par
    current <- accepted$fit
    r <- accepted$residual
    remember(theta,current)
    last.search <- accepted$search
    direction.counts[accepted$kind] <- direction.counts[accepted$kind] + 1L
  }
  converged <- max(abs(root.residual(current))) <= thres
  if (converged) reason <- "adjusted-score tolerance reached"
  if(!attempt.recorded) record.attempt(reason)
  ## Failed attempts remain failed. Archive the best finite candidate found,
  ## without turning its residual minimum into a converged estimate.
  if(!converged) {
    theta <- best$par; current <- best$fit; working.scale <- best$working.scale
  }
  list(par = theta, convergence = converged, step = step,
       diagnostics = list(method = "fixed_weight_mean_br", solver = "joint_damped_newton_armijo_svd",
                          reason = reason, iterations = step, evaluations = evaluations,
                          max_scaled_score = max(abs(root.residual(current))), score_tolerance = thres,
                          max_abs_score = max(abs(current$adjusted.score)),
                          adjusted_score = current$adjusted.score, score_scale = score.scale,
                          working_score_scale=working.scale,
                          weights_normalized_to_mean_one = TRUE,
                          direction_steps=direction.counts, line_searches=line.searches,
                          rejected_steps=rejected.steps, invalid_candidates=invalid.candidates,
                          maximum_jacobian_condition=max.condition,
                          last_residual_stationarity=last.stationarity,
                          last_line_search=last.search, last_candidate_error=last.error,
                          restarts=restarts, restart_reason=restart.reason,
                          restart_history=restart.history,
                          additional_start_max_iterations=50L,
                          maximum_iterations=max.step))
}
