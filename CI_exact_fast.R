exact <- function(param, y, x, va, vb, weight=NULL,
                  max.step, thres=1e-3, thres.dicho=1e-3,
                  pars, se, pa, pb, optim.maxit = 50,
                  optim.reltol = 1e-6){
  
  if (is.null(weight)) {
    weight <- rep(1,length(y))
  }
  
  ## ------------------------------------------------------------
  ## Setup
  ## ------------------------------------------------------------
  getProb <- if (param == "RR") getProbRR else getProbRD
  
  alpha.ml <- pars[1:pa]
  beta.ml  <- pars[(pa + 1):(pa + pb)]
  
  ## Precompute indices to avoid repeated x == 0 / x == 1
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
  n  <- length(y)
  
  eps <- 1e-12
  
  ## Warm start: use ML estimates instead of zeros
  alpha.start <- alpha.ml
  beta.start  <- beta.ml
  
  ## ------------------------------------------------------------
  ## Negative log-likelihood helper
  ## ------------------------------------------------------------
  nll_fun <- function(alpha, beta, y.local){
    
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
    lower <- rep(-bound, length(par))
    upper <- rep(bound, length(par))
    
    fit <- tryCatch(
      stats::optim(
        par,
        fn,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = optim.maxit, factr = optim.reltol / .Machine$double.eps)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit) || any(!is.finite(fit$par))) {
      return(list(par = par, value = fn(par), convergence = 99))
    }
    
    fit
  }
  
  ## ------------------------------------------------------------
  ## Profile nuisance optimization
  ## ------------------------------------------------------------
  optm.beta <- function(alphaj, j, y.local){
    
    alpha <- alpha.start
    beta  <- beta.start
    
    alpha[j] <- alphaj
    
    diff <- thres + 1
    step <- 0
    
    neg.log.likelihood.alpha <- function(alpha.in){
      alpha.in[j] <- alphaj
      nll_fun(alpha.in, beta, y.local)
    }
    
    neg.log.likelihood.beta <- function(beta.in){
      nll_fun(alpha, beta.in, y.local)
    }
    
    while(diff > thres && step < max.step){
      
      step <- step + 1
      
      opt1 <- safe_optim(alpha, neg.log.likelihood.alpha,8)
      
      diff1 <- Diff(opt1$par, alpha)
      alpha <- opt1$par
      alpha[j] <- alphaj
      
      opt2 <- safe_optim(beta, neg.log.likelihood.beta,10)
      
      diff2 <- Diff(opt2$par, beta)
      beta <- opt2$par
      
      diff <- max(diff1, diff2)
    }
    
    nll_fun(alpha, beta, y.local)
  }
  
  ## Cached LRT
  LRT.alpha <- function(alphaj, j, y.local, ll.null = NULL){
    
    if (is.null(ll.null)) {
      ll.null <- optm.beta(alpha.ml[j], j, y.local)
    }
    
    ll.alt <- optm.beta(alphaj, j, y.local)
    
    2 * (ll.alt - ll.null)
  }
  
  ## ------------------------------------------------------------
  ## Simulate distribution of profile-LRT statistic
  ## ------------------------------------------------------------
  ptail <- function(alphaj, j, nsim = 500){
    
    ## First fit nuisance parameter under alpha_j using observed y
    alpha.sim <- alpha.start
    beta.sim  <- beta.start
    
    alpha.sim[j] <- alphaj
    
    diff <- thres + 1
    step <- 0
    
    neg.log.likelihood.alpha.sim <- function(alpha.in){
      alpha.in[j] <- alphaj
      nll_fun(alpha.in, beta.sim, y)
    }
    
    neg.log.likelihood.beta.sim <- function(beta.in){
      nll_fun(alpha.sim, beta.in, y)
    }
    
    while(diff > thres && step < max.step){
      
      step <- step + 1
      
      opt1 <- safe_optim(alpha.sim, neg.log.likelihood.alpha.sim,8)
      
      diff1 <- Diff(opt1$par, alpha.sim)
      alpha.sim <- opt1$par
      alpha.sim[j] <- alphaj
      
      opt2 <- safe_optim(beta.sim, neg.log.likelihood.beta.sim,10)
      
      diff2 <- Diff(opt2$par, beta.sim)
      beta.sim <- opt2$par
      
      diff <- max(diff1, diff2)
    }
    
    ## Fitted probabilities under constrained alpha_j
    prob <- getProb(
      mat_vec_mul(va, alpha.sim),
      mat_vec_mul(vb, beta.sim)
    )
    
    p0 <- prob[idx0, 1]
    p1 <- prob[idx1, 2]
    
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    p1 <- pmin(pmax(p1, eps), 1 - eps)
    
    LRT.sim <- numeric(nsim)
    
    for(i in seq_len(nsim)){
      
      y.sim <- numeric(n)
      
      y.sim[idx0] <- rbinom(n0, 1, p0)
      y.sim[idx1] <- rbinom(n1, 1, p1)
      
      ## Cache null likelihood for this simulated dataset
      ll.null.sim <- optm.beta(alpha.ml[j], j, y.sim)
      
      LRT.sim[i] <- LRT.alpha(
        alphaj,
        j,
        y.sim,
        ll.null = ll.null.sim
      )
    }
    
    LRT.sim
  }
  
  ## ------------------------------------------------------------
  ## Faster acceptability function
  ## ------------------------------------------------------------
  acceptability <- function(alphaj, LRT.obs, LRT.sim){
    
    p.left.obs  <- mean(LRT.sim <= LRT.obs)
    p.right.obs <- mean(LRT.sim >= LRT.obs)
    p.min.obs   <- min(p.left.obs, p.right.obs)
    
    p.left  <- vapply(LRT.sim, function(z) mean(LRT.sim <= z), numeric(1))
    p.right <- vapply(LRT.sim, function(z) mean(LRT.sim >= z), numeric(1))
    p.min   <- pmin(p.left, p.right)
    
    mean(p.min <= p.min.obs)
  }
  
  ## ------------------------------------------------------------
  ## Dichotomy
  ## ------------------------------------------------------------
  dichotomy <- function(j, alpha.low, alpha.up,
                        direction = "low",
                        thres.dicho = 1e-3,
                        max.step = 20,
                        grid.n = 25,
                        nsim.grid = 1500,
                        nsim.refine = 1500){
    
    eval_accept <- function(a, nsim){
      LRT.obs <- LRT.alpha(a, j, y)
      LRT.sim <- ptail(a, j, nsim = nsim)
      acceptability(a, LRT.obs, LRT.sim)
    }
    
    ## Scan outward from alpha.ml[j]
    grid <- if(direction == "low"){
      seq(alpha.up, alpha.low, length.out = grid.n)
    } else {
      seq(alpha.low, alpha.up, length.out = grid.n)
    }
    
    a.vals <- vapply(grid, eval_accept, numeric(1), nsim = nsim.grid)
    accepted <- a.vals >= 0.05
    
    ## If even the ML-side point is rejected, something is numerically wrong
    if(!accepted[1]){
      return(list(
        alpha.dicho = grid[1],
        convergence = FALSE,
        reason = "MLE-side point rejected",
        grid = grid,
        a.vals = a.vals
      ))
    }
    
    ## Find the first pair of consecutive rejections moving outward.
    ## This is more stable than using the first rejection when the
    ## acceptability curve is noisy in small rare-event samples.
    rejected <- !accepted
    consecutive.reject <- which(rejected[-length(rejected)] & rejected[-1]) + 1
    
    ## If no consecutive rejection found, CI endpoint is beyond search range
    if(length(consecutive.reject) == 0){
      return(list(
        alpha.dicho = grid[length(grid)],
        convergence = FALSE,
        reason = "consecutive rejected points not found in search range",
        grid = grid,
        a.vals = a.vals
      ))
    }
    
    outer.idx <- consecutive.reject[1]
    inner.idx <- max(which(accepted[seq_len(outer.idx)]))
    
    ## Bracket: inner accepted, outer rejected
    inner <- grid[inner.idx]
    outer <- grid[outer.idx]
    
    step <- 0
    
    while(abs(inner - outer) > thres.dicho && step < max.step){
      
      step <- step + 1
      mid <- (inner + outer) / 2
      
      a.mid <- eval_accept(mid, nsim = nsim.refine)
      
      if(a.mid >= 0.05){
        inner <- mid
      } else {
        outer <- mid
      }
    }
    
    list(
      alpha.dicho = inner,
      convergence = step < max.step,
      grid = grid,
      a.vals = a.vals
    )
  }
  
  ## ------------------------------------------------------------
  ## Get candidate alpha bounds
  ## -----------------------------------------------------------

  cilen <- ifelse(length(y) < 150, 8, 4)
  alpha.up.start  <- pmin(alpha.ml + cilen * se[1:pa], 8)
  alpha.low.start <- pmax(alpha.ml - cilen * se[1:pa], -8)
  #alpha.up.start  <- 8
  #alpha.low.start <- -8
  
  alpha.up1 <- rep(0, pa)
  
  for(j in seq_len(pa)){
    alpha.up1[j] <- dichotomy(
      j,
      alpha.ml[j],
      alpha.up.start[j],
      direction = "up",
      thres.dicho = thres.dicho
    )$alpha.dicho
  }
  
  alpha.low1 <- rep(0, pa)
  
  for(j in seq_len(pa)){
    alpha.low1[j] <- dichotomy(
      j,
      alpha.low.start[j],
      alpha.ml[j],
      direction = "low",
      thres.dicho = thres.dicho
    )$alpha.dicho
  }
  
  ## ------------------------------------------------------------
  ## P-values
  ## ------------------------------------------------------------
  p.value <- rep(0, pa)
  
  for(j in seq_len(pa)){
    
    LRT.obs.p <- LRT.alpha(0, j, y)
    
    LRT.sim.p <- ptail(
      0,
      j,
      nsim = 1500                       ##### number of simulation
    )
    
    p.value[j] <- acceptability(
      0,
      LRT.obs.p,
      LRT.sim.p
    )
  }
  
  list(
    low = alpha.low1,
    up  = alpha.up1,
    p   = p.value
  )
}
