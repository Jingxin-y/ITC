#' Maximum Likelihood Estimation via Alternating Optimization

#' model by alternately optimizing the negative log-likelihood over \code{alpha} and \code{beta} until convergence or until a maximum number of iterations is reached.
#'
#' @param param Character. Specifies the model type: use \code{"RR"} for relative risk , use \code{"RD"} for risk difference.
#' @param y Numeric vector. Binary response variable (0/1) of length \code{n}.
#' @param x Numeric. Group indicator (0/1) of length \code{n}.
#' @param va Numeric matrix. Design matrix for the \code{alpha} parameters (dimensions \code{n} × \code{pa}).
#' @param vb Numeric matrix. Design matrix for the \code{beta} parameters (dimensions \code{n} × \code{pb}).
#' @param alpha.start Numeric vector of length \code{pa}. Initial values for the \code{alpha} parameters.
#' @param beta.start Numeric vector of length \code{pb}. Initial values for the \code{beta} parameters.
#' @param weight Numeric vector of length \code{n}. Observation weight.
#' @param max.step Integer. Maximum number of alternating optimization iterations.
#' @param thres Numeric. Convergence threshold: the algorithm stops when the relative change in parameters falls below this value.
#' @param pa Integer. Number of \code{alpha} parameters (length of \code{alpha.start}).
#' @param pb Integer. Number of \code{beta} parameters (length of \code{beta.start}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{par}}{Numeric vector of length \code{pa + pb}: the estimated parameters \code{c(alpha, beta)}.}
#'   \item{\code{convergence}}{Logical. \code{TRUE} if the algorithm converged within \code{max.step} iterations; otherwise \code{FALSE}.}
#'   \item{\code{value}}{Numeric. Negative log-likelihood evaluated at the final parameter estimates.}
#'   \item{\code{step}}{Integer. Number of iterations actually performed.}
#' }
#'
max.likelihood = function(param, y, x, va, vb, alpha.start, beta.start, weight, 
                          max.step, thres, pa, pb) {
    
    startpars = c(alpha.start, beta.start)
    
    getProb = if (param == "RR") getProbRR else getProbRD
    
    ## negative log likelihood function
    neg.log.likelihood = function(pars) {
        alpha = pars[1:pa]
        beta = pars[(pa + 1):(pa + pb)]
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0 = p0p1[, 1];   p1 = p0p1[, 2]
        
        return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] + 
                        (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x == 
                                                                                          1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x == 
                                                                                                                                                                 1]) * weight[x == 1]))
    }
    
    neg.log.likelihood.alpha = function(alpha){
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                        (y[x==0])*log(p0[x==0])*weight[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                           (y[x==1])*log(p1[x==1])*weight[x==1]))  
    }
    
    neg.log.likelihood.beta = function(beta){
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                        (y[x==0])*log(p0[x==0])*weight[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                           (y[x==1])*log(p1[x==1])*weight[x==1]))  
    }
    
    
    ## Optimization 
    
    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    alpha = alpha.start; beta = beta.start
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
        step = step + 1
        opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
        diff1 = Diff(opt1$par,alpha)
        alpha = opt1$par
        opt2 = stats::optim(beta,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
        diff  = max(diff1,Diff(opt2$par,beta))
        beta = opt2$par
    }
    
    opt = list(par = c(alpha,beta), convergence = (step < max.step), 
               value = neg.log.likelihood(c(alpha,beta)), step = step)
    
    return(opt)
}

