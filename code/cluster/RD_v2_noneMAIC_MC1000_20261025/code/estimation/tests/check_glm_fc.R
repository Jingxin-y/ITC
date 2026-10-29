## Independent closed-form and numerical Taylor checks; run from anchored_v1.
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
close_to <- function(a, b, tol = 1e-6) stopifnot(isTRUE(all.equal(a, b,
  tolerance = tol, check.attributes = FALSE)))

## Two saturated arms, including zero events. Nonuniform analytic weights
## have correction k/2, where k=sum(w^2)/sum(w), rather than 1/2 replicated counts.
dat <- data.frame(x = rep(0:1, each = 20),
  y = c(rep(0, 20), rep(1, 5), rep(0, 15)))
for (w in list(rep(1, 40), rep(c(.3, .7, 1.5, 2), 10))) {
  wn <- w/mean(w)
  for (spec in list(c("binomial", "bernoulli"), c("poisson", "poisson"),
                    c("poisson", "bernoulli"))) {
    z <- fit_logglm_fc(dat, w, family = spec[1], bias_model = spec[2],
      model_covariates = character(0), tolerance = 1e-9)
    stopifnot(z$converged)
    expected <- vapply(0:1, function(arm) {
      ix <- dat$x == arm; W <- sum(wn[ix]); S <- sum(wn[ix]*dat$y[ix])
      k <- sum(wn[ix]^2)/W
      (S + k/2)/(W + if (spec[2] == "bernoulli") k/2 else 0)
    }, numeric(1))
    close_to(predict(z$fit, data.frame(x = 0:1), type = "response"), expected, 1e-7)
    scaled <- fit_logglm_fc(dat, w*17, family = spec[1], bias_model = spec[2],
      model_covariates = character(0), tolerance = 1e-9)
    close_to(z$point.est, scaled$point.est, 1e-8); close_to(z$vcov, scaled$vcov, 1e-8)
    X <- model.matrix(z$fit); mu <- fitted(z$fit)
    h <- if (spec[1] == "binomial") mu/(1-mu) else mu
    v <- if (spec[1] == "binomial") mu/(1-mu) else
      if (spec[2] == "poisson") mu else mu*(1-mu)
    Ji <- solve(crossprod(X, X*(wn*h)))
    close_to(z$vcov, Ji %*% crossprod(X, X*(wn^2*v)) %*% Ji)
  }
}
cat("PASS: independent saturated-arm solutions, zero events, covariance and weight scaling\n")

## Numerical random-Taylor oracle: differentiate each unadjusted subject score.
## For Poisson, a two-node moment rule reproduces E(Y) and Var(Y); these score
## derivatives are affine in Y, so the required moments are exact.
set.seed(1301)
X <- cbind(1, rep(0:1, 6), runif(12, -.8, .8)); theta <- c(-2, .25, -.2)
w <- runif(12, .3, 2); mu <- exp(drop(X %*% theta)); p <- ncol(X)
oracle <- function(family, bias_model) {
  h <- 2e-4; J <- B <- matrix(0, p, p); second <- array(0, c(p, p, p)); terms <- list()
  for (i in seq_len(nrow(X))) {
    values <- if (bias_model == "bernoulli") 0:1 else mu[i]+c(-1,1)*sqrt(mu[i])
    probs <- if (bias_model == "bernoulli") c(1-mu[i],mu[i]) else c(.5,.5)
    for (k in 1:2) {
      u <- function(b) {
        m <- exp(sum(X[i,]*b)); r <- values[k]-m
        w[i]*X[i,]*(if (family == "binomial") r/(1-m) else r)
      }
      u0 <- u(theta); first <- matrix(0,p,p); d2 <- array(0,c(p,p,p))
      for (j in seq_len(p)) {
        ej <- numeric(p); ej[j] <- h
        first[,j] <- (u(theta+ej)-u(theta-ej))/(2*h)
        for (l in seq_len(p)) {
          el <- numeric(p); el[l] <- h
          d2[,j,l] <- if (j == l) (u(theta+ej)-2*u0+u(theta-ej))/h^2 else
            (u(theta+ej+el)-u(theta+ej-el)-u(theta-ej+el)+u(theta-ej-el))/(4*h^2)
        }
      }
      J <- J-probs[k]*first; B <- B+probs[k]*tcrossprod(u0)
      second <- second+probs[k]*d2
      terms[[length(terms)+1L]] <- list(prob=probs[k],first=first,score=u0)
    }
  }
  Ji <- solve(J); V <- Ji %*% B %*% t(Ji)
  product <- Reduce(`+`,lapply(terms,function(z) z$prob*drop(z$first %*% Ji %*% z$score)))
  contraction <- vapply(seq_len(p),function(j) sum(second[j,,]*V),numeric(1))
  -product-.5*contraction
}
for (spec in list(c("binomial","bernoulli"),c("poisson","poisson"),c("poisson","bernoulli"))) {
  z <- glm_fc_score(theta,X,rep(0,12),w,spec[1],spec[2])
  close_to(z$adjustment,oracle(spec[1],spec[2]),5e-6)
}
cat("PASS: all three adjustments agree with an independent numerical Taylor oracle\n")

## Unit-weight Poisson adjustment is the gradient of half log(det(I)).
penalty <- function(b) .5*as.numeric(determinant(crossprod(X,X*exp(drop(X%*%b))),logarithm=TRUE)$modulus)
grad <- vapply(seq_len(p),function(j) { e <- numeric(p); e[j] <- 1e-5
  (penalty(theta+e)-penalty(theta-e))/2e-5 },numeric(1))
close_to(glm_fc_score(theta,X,rep(0,12),rep(1,12),"poisson","poisson")$adjustment,grad,1e-7)
lb <- glm_fc_score(theta,X,rep(0,12),rep(1,12),"binomial","bernoulli")
hat <- mu/(1-mu)*rowSums((X%*%solve(lb$information))*X)
close_to(lb$adjustment,drop(crossprod(X,.5*hat)),1e-8)
cat("PASS: unit weights recover standard Poisson-Firth and log-binomial AS_mean equations\n")

## Robust covariance independently assembled from the unadjusted score.
dat$z <- rep(c(-.5,.5),20)
dat$y <- as.integer(seq_len(40) %% 5 == 0)
z <- fit_robust_lp_fc(dat,model_covariates="z")
stopifnot(z$converged)
XX <- model.matrix(z$fit); pp <- fitted(z$fit); Ji <- solve(crossprod(XX,XX*pp))
expected <- Ji %*% crossprod(XX*(dat$y-pp)) %*% Ji
close_to(z$vcov,expected,1e-8)

## Singular designs are rejected; exhausted iterations are not valid estimates.
bad <- tryCatch(fit_lp_fc(transform(dat,z=x),model_covariates="z"),error=function(e)e)
stopifnot(inherits(bad,"error"))
failed <- suppressWarnings(fit_lp_fc(dat,model_covariates="z",maxit=1,tolerance=1e-14))
stopifnot(!failed$converged,all(is.na(failed$vcov)))
rows <- comparison_parameter_rows("LP-FC",failed,"robust_glm")
stopifnot(all(rows$status=="fit_not_converged"),all(is.na(rows$estimate)))
success <- transform(dat,y=1)
boundary <- suppressWarnings(fit_logbinomial_fc(success,model_covariates=character()))
stopifnot(!boundary$converged,all(is.na(boundary$vcov)))
cat("PASS: empirical covariance, singular design, iteration failure and boundary failure\n")
