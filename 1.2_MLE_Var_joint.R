
#' Calculate the variance–covariance matrix of MLEs for the relative risk (RR) model and the riskdifference (RD) model
#'
#' This function computes the inverse of the observed Fisher information matrix

#'
#' @param x numeric vector of length \code{n}. Binary exposure indicator (0/1).
#' @param alpha.ml Numeric vector of length \code{p_a}. Estimated \code{\alpha} parameters.
#' @param beta.ml Numeric vector of length \code{p_b}. Estimated \code{\beta} parameters.
#' @param va Numeric matrix of dimension \code{n \times p_a}. Design matrix for the \code{\alpha} component.
#' @param vb Numeric matrix of dimension \code{n \times p_b}. Design matrix for the \code{\beta} component.
#' @param weight Numeric vector of length \code{n}. Observation weight.
#'
#' @return A \code{p}-by-\code{p} matrix (\code{p = length(alpha.ml) + length(beta.ml)}), the variance–covariance matrix of the MLEs).
#'
#' need the package "MASS"


make_target_mean_vcov <- function(target_var = NULL,
                                  m_AD = NULL,
                                  target_cov = NULL,
                                  target_mean_vcov = NULL) {
  if (!is.null(target_mean_vcov)) {
    return(as.matrix(target_mean_vcov))
  }

  if (is.null(m_AD)) {
    stop("m_AD must be supplied if target_mean_vcov is not supplied.")
  }

  if (!is.null(target_cov)) {
    target_cov <- as.matrix(target_cov)
    return(target_cov / m_AD)
  }

  if (is.null(target_var)) {
    stop("Either target_var, target_cov, or target_mean_vcov must be supplied.")
  }

  if (is.matrix(target_var)) {
    return(as.matrix(target_var) / m_AD)
  }

  target_var <- as.numeric(target_var)
  diag(target_var / m_AD, length(target_var))
}

joint_sandwich_summaryAD <- function(D,
                                     weight,
                                     Hgamma,
                                     Si,
                                     hessian,
                                     u1i,
                                     u2i,
                                     target_var = NULL,
                                     m_AD = NULL,
                                     target_cov = NULL,
                                     target_mean_vcov = NULL,
                                     return_full = FALSE) {
  q <- ncol(D)
  p <- ncol(Si)

  V_mu <- make_target_mean_vcov(
    target_var = target_var,
    m_AD = m_AD,
    target_cov = target_cov,
    target_mean_vcov = target_mean_vcov
  )

  if (!all(dim(V_mu) == c(q, q))) {
    stop("Dimension mismatch: Var(target_mean_hat) must be q by q, where q = ncol(vg).")
  }

  A_mu_mu <- diag(q)
  A_mu_gamma <- matrix(0, q, q)
  A_mu_theta <- matrix(0, q, p)

  A_gamma_mu <- sum(weight) * diag(q)
  A_gamma_gamma <- -t(D) %*% (weight * Hgamma)
  A_gamma_theta <- matrix(0, q, p)

  A_theta_mu <- matrix(0, p, q)
  A_theta_gamma <- -t(Si) %*% (weight * Hgamma)
  A_theta_theta <- hessian

  A <- rbind(
    cbind(A_mu_mu,     A_mu_gamma,     A_mu_theta),
    cbind(A_gamma_mu,  A_gamma_gamma,  A_gamma_theta),
    cbind(A_theta_mu,  A_theta_gamma,  A_theta_theta)
  )

  Ui_ipd <- cbind(u1i, u2i)
  B_ipd <- t(Ui_ipd) %*% Ui_ipd

  B <- matrix(0, nrow = q + q + p, ncol = q + q + p)

  idx_mu <- seq_len(q)
  idx_ipd <- (q + 1):(q + q + p)

  B[idx_mu, idx_mu] <- V_mu
  B[idx_ipd, idx_ipd] <- B_ipd

  Ainv <- MASS::ginv(A)
  V <- Ainv %*% B %*% t(Ainv)
  V <- (V + t(V)) / 2

  idx_gamma <- (q + 1):(2 * q)
  idx_theta <- (2 * q + 1):(2 * q + p)

  if (return_full) {
    return(list(
      full = V,
      target_mean = V[idx_mu, idx_mu, drop = FALSE],
      gamma = V[idx_gamma, idx_gamma, drop = FALSE],
      theta = V[idx_theta, idx_theta, drop = FALSE]
    ))
  }

  V[idx_theta, idx_theta, drop = FALSE]
}

has_target_mean_vcov <- function(target_var = NULL,
                                 m_AD = NULL,
                                 target_cov = NULL,
                                 target_mean_vcov = NULL) {
  !is.null(target_mean_vcov) ||
    (!is.null(m_AD) && (!is.null(target_cov) || !is.null(target_var)))
}

### variance calculation

var.mle.rr.joint = function(y, x, alpha.ml, beta.ml, gamma.ml,
                      va, vb, vg, target_mean,
                      target_var = NULL,
                      m_AD = NULL,
                      target_cov = NULL,
                      target_mean_vcov = NULL) {

  D <- sweep(vg, 2, target_mean, "-")
  weight <- normalize_exp_weights(D %*% gamma.ml)
  Dbar.w <- colSums(weight * D) / sum(weight)
  Hgamma <- sweep(D, 2, Dbar.w, "-")
  
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
    
    
    ### first 
    dpsi0.by.dphi = (1 - p0) * (1 - p1)/((1 - p0) + (1 - p1))
    dpsi0.by.dtheta = -(1 - p0)/((1 - p0) + (1 - p1))
    tmp = cbind((dpsi0.by.dtheta + x) * va, dpsi0.by.dphi * vb)
    ## since dtheta.by.dalpha = va, and dphi.by.dbeta = vb
    ### second
    
    ### Building blocks
    
    dpsi0_by_dtheta <- -(1 - p0) / (1 - p0 + 1 - p1)
    dpsi0_by_dphi <- (1 - p0) * (1 - p1) / (1 - p0 + 1 - p1)
    
    dtheta_by_dalpha <- va
    dphi_by_dbeta <- vb
    
    dl_by_dpsi0 <- (y - pA) / (1 - pA)
    d2l_by_dpsi0_2 <- (y - 1) * pA / ((1 - pA)^2)
    
    
    ###### d2l_by_dalpha_2
    
    d2psi0_by_dtheta_2 <- ((p0 - p1) * dpsi0_by_dtheta - (1 - p0) * p1) / ((1 -
                                                                              p0 + 1 - p1)^2)
    
    d2l_by_dtheta_2 <- d2l_by_dpsi0_2 * (dpsi0_by_dtheta + x)^2 + dl_by_dpsi0 *
      d2psi0_by_dtheta_2
    
    d2l_by_dalpha_2 <- t(dtheta_by_dalpha * d2l_by_dtheta_2 * weight) %*%
      dtheta_by_dalpha
    
    
    ###### d2l_by_dalpha_dbeta
    
    d2psi0_by_dtheta_dphi <- (1 - p0) * (1 - p1) * (p0 - p1) / (1 - p0 + 1 -
                                                                  p1)^3
    
    d2l_by_dtheta_dphi <- d2l_by_dpsi0_2 * (dpsi0_by_dtheta + x) * dpsi0_by_dphi +
      dl_by_dpsi0 * d2psi0_by_dtheta_dphi
    
    d2l_by_dalpha_dbeta <- t(dtheta_by_dalpha * d2l_by_dtheta_dphi * weight) %*%
      dphi_by_dbeta
    d2l_by_dbeta_dalpha <- t(d2l_by_dalpha_dbeta)
    # d2l_by_dalpha_dbeta is symmetric itself if (because) va=vb
    
    
    #### d2l_by_dbeta2
    
    d2psi0_by_dphi_2 <- (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2) / (1 - p0 + 1 -
                                                                   p1)^2) * dpsi0_by_dphi
    
    d2l_by_dphi_2 <- d2l_by_dpsi0_2 * (dpsi0_by_dphi)^2 + dl_by_dpsi0 * d2psi0_by_dphi_2
    
    d2l_by_dbeta_2 <- t(dphi_by_dbeta * d2l_by_dphi_2 * weight) %*% dphi_by_dbeta
    
    
    hessian <- -rbind(cbind(d2l_by_dalpha_2, d2l_by_dalpha_dbeta), cbind(
      d2l_by_dbeta_dalpha,
      d2l_by_dbeta_2
    ))
    
    ### blocks
    u1i <- weight * D
    dl_by_dpsi0 <- (y - pA) / (1 - pA)
    u2i <- (weight * dl_by_dpsi0) * tmp
    
    
    ### martrix A
    A11 <- -t(D) %*% (weight * Hgamma)
    A12 <- matrix(0, ncol(D), ncol(tmp))
    Si <- dl_by_dpsi0 * tmp
    A21 <- -t(Si) %*% (weight * Hgamma)
    A22 <- hessian
    A <- rbind(
      cbind(A11, A12),
      cbind(A21, A22)
    )
    
    if (has_target_mean_vcov(target_var, m_AD, target_cov, target_mean_vcov)) {
      V.full <- joint_sandwich_summaryAD(
        D = D,
        weight = weight,
        Hgamma = Hgamma,
        Si = Si,
        hessian = hessian,
        u1i = u1i,
        u2i = u2i,
        target_var = target_var,
        m_AD = m_AD,
        target_cov = target_cov,
        target_mean_vcov = target_mean_vcov,
        return_full = TRUE
      )$full

      q <- ncol(D)
      idx <- (q + 1):ncol(V.full)
      return(V.full[idx, idx, drop = FALSE])
    }

    ### matrix B
    Ui <- cbind(u1i, u2i)
    B <- t(Ui) %*% Ui
    V <- ginv(A) %*% B %*% t(ginv(A))

    return(V)
}




### variance calculation

var.mle.rd.joint = function(y, x, alpha.ml, beta.ml, gamma.ml,
                      va, vb, vg, target_mean,
                      target_var = NULL,
                      m_AD = NULL,
                      target_cov = NULL,
                      target_mean_vcov = NULL) {
  D <- sweep(vg, 2, target_mean, "-")
  weight <- normalize_exp_weights(D %*% gamma.ml)
  Dbar.w <- colSums(weight * D) / sum(weight)
  Hgamma <- sweep(D, 2, Dbar.w, "-")

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

    dp0.by.dphi = s0 * s1/(s0 + s1)
    dp0.by.drho = -s0/(s0 + s1)
    drho.by.dalpha = (1 - rho^2) * va
    dphi.by.dbeta = vb

    tmp = cbind((dp0.by.drho + x) * drho.by.dalpha, dp0.by.dphi * dphi.by.dbeta)
    

    ### First order derivatives ###
    
    dl_by_dpA <- (y - pA) / sA
    dp0_by_dphi <- s0 * s1 / (s0 + s1)
    dp0_by_drho <- -s0 / (s0 + s1)
    drho_by_dalpha <- va * (1 - rho^2)
    dphi_by_dbeta <- vb
    
    dpA_by_drho <- dp0_by_drho + x
    dpA_by_dalpha <- drho_by_dalpha * dpA_by_drho
    dpA_by_dphi <- dp0_by_dphi
    dpA_by_dbeta <- dphi_by_dbeta * dpA_by_dphi
    
    ### Second order derivatives ###
    
    d2l_by_dpA_2 <- -(y - pA)^2 / sA^2
    d2pA_by_drho_2 <- s0 * s1 * (2 - 2 * p0 - 2 * p1) / (s0 + s1)^3
    d2pA_by_dphi_drho <- (s0 * (1 - 2 * p1) - s1 * (1 - 2 * p0)) * s0 * s1 / (s0 +
                                                                                s1)^3
    d2pA_by_dphi_2 <- (s0^2 * (1 - 2 * p1) + s1^2 * (1 - 2 * p0)) * s0 * s1 / (s0 +
                                                                                 s1)^3
    
    d2rho_by_dalpha_2 <- -2 * t(va * rho) %*% drho_by_dalpha
    
    ### Compute elements of the Hessian matrix ###
    
    d2l_by_dalpha_2 <- t(dpA_by_dalpha * d2l_by_dpA_2 * weight) %*% dpA_by_dalpha +
      t(drho_by_dalpha * dl_by_dpA * d2pA_by_drho_2 * weight) %*% drho_by_dalpha -
      2 * t(va * rho * dl_by_dpA * dpA_by_drho * weight) %*% drho_by_dalpha
    
    d2l_by_dalpha_dbeta <- t(dpA_by_dalpha * d2l_by_dpA_2 * weight) %*% dpA_by_dbeta +
      t(drho_by_dalpha * dl_by_dpA * d2pA_by_dphi_drho * weight) %*% dphi_by_dbeta
    d2l_by_dbeta_dalpha <- t(d2l_by_dalpha_dbeta)
    
    d2l_by_dbeta_2 <- t(dpA_by_dbeta * d2l_by_dpA_2 * weight) %*% dpA_by_dbeta +
      t(dphi_by_dbeta * dl_by_dpA * d2pA_by_dphi_2 * weight) %*% dphi_by_dbeta
    
    hessian <- -rbind(cbind(d2l_by_dalpha_2, d2l_by_dalpha_dbeta), cbind(
      d2l_by_dbeta_dalpha,
      d2l_by_dbeta_2
    ))
    
    u1i <- weight * D
    u2i <- (weight * dl_by_dpA) * tmp
    
    
    ### martrix A
    A11 <- -t(D) %*% (weight * Hgamma)
    A12 <- matrix(0, ncol(D), ncol(tmp))
    Si <- dl_by_dpA * tmp
    A21 <- -t(Si) %*% (weight * Hgamma)
    A22 <- hessian
    A <- rbind(
      cbind(A11, A12),
      cbind(A21, A22)
    )

    if (has_target_mean_vcov(target_var, m_AD, target_cov, target_mean_vcov)) {
      V.full <- joint_sandwich_summaryAD(
        D = D,
        weight = weight,
        Hgamma = Hgamma,
        Si = Si,
        hessian = hessian,
        u1i = u1i,
        u2i = u2i,
        target_var = target_var,
        m_AD = m_AD,
        target_cov = target_cov,
        target_mean_vcov = target_mean_vcov,
        return_full = TRUE
      )$full

      q <- ncol(D)
      idx <- (q + 1):ncol(V.full)
      return(V.full[idx, idx, drop = FALSE])
    }

    ### matrix B
    Ui <- cbind(u1i, u2i)
    B <- t(Ui) %*% Ui
    V <- ginv(A) %*% B %*% t(ginv(A))
    
    return(V)
}


