#' Bayesian estimator for Risk Difference on Fisher-z scale
#'
#' @description
#' Conjugate-Beta model for two independent binomials in treatment (\eqn{x=1}) and control (\eqn{x=0}).
#' Draws \eqn{p_1 \sim \text{Beta}(a_1+N_{1,1},\, b_1+N_{a1}-N_{1,1})} and
#' \eqn{p_0 \sim \text{Beta}(a_0+N_{0,1},\, b_0+N_{a0}-N_{0,1})} by Monte Carlo,
#' forms the risk difference \eqn{d = p_1 - p_0}, and reports summaries on the
#' Fisher-\eqn{z} scale \eqn{\alpha = \operatorname{atanh}(d)} (stabilizes near the
#' boundaries \eqn{[-1,1]}).
#'
#' @param Na0 Integer. Number at risk in control arm (\eqn{x=0}).
#' @param Na1 Integer. Number at risk in treatment arm (\eqn{x=1}).
#' @param N0_1 Integer. Number of events in control arm.
#' @param N1_1 Integer. Number of events in treatment arm.
#' @param a1,b1,a0,b0 Positive numerics. Beta prior hyperparameters for
#'   \eqn{p_1 \sim \text{Beta}(a_1,b_1)} and \eqn{p_0 \sim \text{Beta}(a_0,b_0)}.
#'   Defaults are Jeffreys \code{0.5, 0.5}.
#' @param M Integer. Number of Monte Carlo draws (default \code{1e5}).
#' @param conf Numeric in \code{(0,1)}. Credible mass for intervals (default \code{0.95}).
#'
#' @return A list with components (all on the \eqn{\alpha=\operatorname{atanh}(d)} scale):
#' \describe{
#'   \item{\code{point.est}}{Posterior mean of \eqn{\alpha}.}
#'   \item{\code{se.est}}{Posterior SD of \eqn{\alpha}.}
#'   \item{\code{conf.lower}, \code{conf.upper}}{Equal-tail credible interval endpoints.}
#'   \item{\code{ET}}{Length-2 vector of equal-tail endpoints.}
#'   \item{\code{HPD}}{Length-2 Highest Posterior Density interval from \pkg{HDInterval}.}
#'   \item{\code{p.value}}{Posterior probability \eqn{\Pr(d>0)} (one-sided support for RD>0).}
#' }
#'

bayes_est_RD <- function(Na0, Na1, N0_1, N1_1, a1=.5, b1=.5, a0=.5, b0=.5,
                      M=1e5, conf=0.95){
  p1 <- rbeta(M, a1 + N1_1, b1 + Na1-N1_1)
  p0 <- rbeta(M, a0 + N0_1, b0 + Na0-N0_1)
  d  <- p1 - p0
  alpha <- atanh(d)
  sd <- sd(alpha)
  et <- quantile(alpha, c((1-conf)/2, 1-(1-conf)/2))
  hpd<- HDInterval::hdi(alpha, credMass=conf)
  list(point.est = mean(alpha), se.est = sd, conf.lower = min(et),
       conf.upper = max(et), ET = et, HPD = hpd,
       p.value = mean(d>0))
}

#' Bayesian estimator for Risk Ratio on log scale
#'
#' @description
#' Conjugate-Beta model as above, but summarizes the risk ratio \eqn{d = p_1/p_0}
#' on the log scale \eqn{\alpha=\log d}.
#'
#' @inheritParams bayes_est_RD
#'
#' @return A list with components (all on the \eqn{\alpha=\log(p_1/p_0)} scale):
#' \describe{
#'   \item{\code{point.est}}{Posterior mean of \eqn{\alpha}.}
#'   \item{\code{se.est}}{Posterior SD of \eqn{\alpha}.}
#'   \item{\code{conf.lower}, \code{conf.upper}}{Equal-tail credible interval endpoints.}
#'   \item{\code{ET}}{Equal-tail endpoints.}
#'   \item{\code{HPD}}{HPD interval via \pkg{HDInterval}.}
#'   \item{\code{p.value}}{Posterior probability \eqn{\Pr(d>0)} (for RR>0 which is always true; often
#'         redefine as \eqn{\Pr(\log RR>0)} if desired).}
#' }

bayes_est_RR <- function(Na0, Na1, N0_1, N1_1, a1=.5, b1=.5, a0=.5, b0=.5,
                         M=1e5, conf=0.95){
  p1 <- rbeta(M, a1 + N1_1, b1 + Na1-N1_1)
  p0 <- rbeta(M, a0 + N0_1, b0 + Na0-N0_1)
  d  <- p1/p0
  alpha <- log(d)
  sd <- sd(alpha)
  et <- quantile(alpha, c((1-conf)/2, 1-(1-conf)/2))
  hpd<- HDInterval::hdi(alpha, credMass=conf)
  list(point.est = mean(alpha), se.est = sd, conf.lower = min(et),
       conf.upper = max(et), ET = et, HPD = hpd,
       p.value = mean(d>0))
}



#' g-computation helper functions
#'
#' @description
#' This set of helper functions (\code{mu.est}, \code{l.mu}, \code{var.est},
#' \code{m}, \code{m.prime}, \code{m.prime.prime}, \code{fish}, \code{hii},
#' \code{phi}) are implementations used in the
#' g-computation framework.  
#' The code is written according to the methodology described in
#' \url{https://arxiv.org/pdf/2509.07369}.

mu.est <- function(y,x,beta){mean(c(y,m(x%*%beta)))}


l.mu <- function(y1,x1,beta1,y0,x0,beta0){
    mi <- c(m(x1%*%beta1),m(x0%*%beta0))
    mi1 <- c(m(x1%*%beta1),m(x0%*%beta1))
    mi0 <- c(m(x1%*%beta0),m(x0%*%beta0))
    mu1 <- mu.est(y1,x0,beta1)
    mu0 <- mu.est(y0,x1,beta0)
    
    li1 <- c((1+hii(x1,beta1))*(y1-mi[1:length(y1)]),rep(0,length(y0)))/(length(y1)/(length(c(y1,y0))))+mi1-mu1
    li0 <- c(rep(0,length(y1)),(1+hii(x0,beta0))*(y0-mi[(length(y1)+1):length(mi)]))/(length(y0)/(length(c(y1,y0))))+mi0-mu0
    return(cbind(li0,li1))
  }

var.est <- function(li,p0,p1){
    cov <- var(li)/nrow(li)
    return(cov[1,1]/(p0^2)+cov[2,2]/(p1^2)-2*cov[1,2]/(p0*p1))
  }
  
m <- function(x) plogis(x)
m.prime <- function(x){
  s <- plogis(x)
  s * (1 - s)
}

m.prime.prime <- function(x){
  s <- plogis(x)
  s * (1 - s) * (1 - 2*s)
}


fish <- function(x, beta){
  eta <- drop(x %*% beta)
  w <- m.prime(eta)                 # length n
  # t(x) %*% diag(w) %*% x  等价于  crossprod(x * w, x)
  crossprod(x * w, x) / nrow(x)
}
hii <- function(x,beta){m.prime(x%*%beta)*rowSums(x%*%ginv(nrow(x)*fish(x,beta))*x)}
phi <- function(y,x,beta,p){x%*%ginv(fish(x,beta))*as.vector(y-m(x%*%beta))/p}

#' Transform point/SE/CI for RD to Fisher-z scale
#'
#' @description
#' Converts an estimate on the RD scale (\eqn{d}) and its SE/CI to the stabilized
#' Fisher-\eqn{z} scale \eqn{\alpha=\operatorname{atanh}(d)} using the delta method:
#' \eqn{\mathrm{SE}(\alpha) \approx \mathrm{SE}(d)/(1-d^2)}; CI endpoints are
#' transformed via \code{atanh}.
#'
#' @param est Numeric scalar. Point estimate on RD scale.
#' @param se Numeric scalar. Standard error of \code{est} on RD scale.
#' @param conf Numeric length-2 vector. Confidence interval endpoints on RD scale.
#'
#' @return A list with
#' \describe{
#'   \item{\code{point.est}}{\eqn{\operatorname{atanh}(est)}.}
#'   \item{\code{se.est}}{\eqn{se/(1-est^2)}.}
#'   \item{\code{CI}}{Length-2 vector \code{atanh(conf)}.}
#' }
#'
get_estimate <- function(est,se,conf){
  Ealpha <- atanh(est) #+est*se^2/((1-est^2)^2)
  Valpha <- se/(1-est^2)
  CIalpha <- atanh(conf)
  list(point.est = Ealpha, se.est = Valpha, CI = CIalpha)
}

