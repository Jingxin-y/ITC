
#' Calculate the variance–covariance matrix of MLEs for the relative risk (RR) model and the riskdifference (RD) model
#'
#' This function computes the inverse of the observed Fisher information matrix

#'
#' @param x numeric vector of length \code{n}. Binary exposure indicator (0/1).
#' @param alpha.ml Numeric vector of length \code{p_a}. Estimated \code{\alpha} parameters.
#' @param beta.ml Numeric vector of length \code{p_b}. Estimated \code{\beta} parameters.
#' @param va Numeric matrix of dimension \code{n \times p_a}. Design matrix for the \code{\alpha} component.
#' @param vb Numeric matrix of dimension \code{n \times p_b}. Design matrix for the \code{\beta} component.
#' @param weights Numeric vector of length \code{n}. Observation weights.
#'
#' @return A \code{p}-by-\code{p} matrix (\code{p = length(alpha.ml) + length(beta.ml)}), the variance–covariance matrix of the MLEs).
#'
#' need the package "MASS"


### variance calculation

var.mle.rr = function(x, alpha.ml, beta.ml, va, vb, weight) {

    p0p1 = getProbRR(va %*% alpha.ml, vb %*% beta.ml)
    n = dim(vb)[1]
    pA = rep(NA, n)
    p0    = p0p1[,1];  p1 = p0p1[,2]
    pA[x == 0] = p0p1[x == 0, 1]
    pA[x == 1] = p0p1[x == 1, 2]
    eps = 1e-8
    p0 = pmin(pmax(p0, eps), 1 - eps)
    p1 = pmin(pmax(p1, eps), 1 - eps)
    pA = pmin(pmax(pA, eps), 1 - eps)
    p0p1 = cbind(p0, p1)


    expect.dl.by.dpsi0.squared = (pA)/(1 - pA)
    dpsi0.by.dphi = (1 - p0p1[, 1]) * (1 - p0p1[, 2])/((1 - p0p1[, 1]) + (1 -
        p0p1[, 2]))
    dpsi0.by.dtheta = -(1 - p0p1[, 1])/((1 - p0p1[, 1]) + (1 - p0p1[, 2]))
    tmp = cbind((dpsi0.by.dtheta + x) * va, dpsi0.by.dphi * vb)
    ## since dtheta.by.dalpha = va, and dphi.by.dbeta = vb
    fisher.info = (t(expect.dl.by.dpsi0.squared * weight * tmp) %*% tmp)
    return(ginv(fisher.info))
}




### variance calculation

var.mle.rd = function(x, alpha.ml, beta.ml, va, vb, weight) {

    p0p1 = getProbRD(va %*% alpha.ml, vb %*% beta.ml)
    # p0p1 = cbind(p0, p1): n * 2 matrix
    p0 = p0p1[, 1]
    p1 = p0p1[, 2]

    n = nrow(va)
    pA = p0             # P(Y=1|A,V); here A = X
    pA[x == 1] = p1[x == 1]
    eps = 1e-8
    p0 = pmin(pmax(p0, eps), 1 - eps)
    p1 = pmin(pmax(p1, eps), 1 - eps)
    pA = pmin(pmax(pA, eps), 1 - eps)
    s0 = p0 * (1 - p0)
    s1 = p1 * (1 - p1)
    sA = pA * (1 - pA)

    rho = as.vector(tanh(va %*% alpha.ml))  #estimated risk differences

    expect.dl.by.dpA.squared = 1/sA
    dp0.by.dphi = s0 * s1/(s0 + s1)
    dp0.by.drho = -s0/(s0 + s1)
    drho.by.dalpha = (1 - rho^2) * va
    dphi.by.dbeta = vb

    tmp = cbind((dp0.by.drho + x) * drho.by.dalpha, dp0.by.dphi * dphi.by.dbeta)
    fisher.info = (t(expect.dl.by.dpA.squared * weight * tmp) %*% tmp)
    return(ginv(fisher.info))
}


