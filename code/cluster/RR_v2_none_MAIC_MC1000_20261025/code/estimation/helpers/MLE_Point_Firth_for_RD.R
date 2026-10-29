## RD: theta=atanh(RD), phi=log(odds product); rho=tanh(theta).
## Derivative source: brmplus commit 3b628a6d32dc7fe0f31fbd3ab181aecd3866f982,
## compare/MLE_Point_Firth_for_RD.R; see BRMPLUS_LICENSE.txt.
## Scalar derivatives for the fixed-weight Bernoulli mean-bias adjustment.
firth_rd_derivatives <- function(alpha, beta, x, va, vb) {
  theta <- as.vector(va %*% alpha)
  probs <- getProbRD(theta, as.vector(vb %*% beta))
  p0 <- pmin(pmax(probs[, 1], 1e-8), 1 - 1e-8)
  p1 <- pmin(pmax(probs[, 2], 1e-8), 1 - 1e-8)
  p <- ifelse(x == 1, p1, p0)
  s0 <- p0 * (1 - p0)
  s1 <- p1 * (1 - p1)
  denom <- s0 + s1
  rho <- tanh(theta)
  rho.first <- 1 - rho^2
  rho.second <- -2 * rho * rho.first

  ## Derivatives of p0 with respect to rho and phi, from the original RD blocks.
  p.rho <- -s0 / denom
  p.phi <- s0 * s1 / denom
  p.rhorho <- s0 * s1 * (2 - 2 * p0 - 2 * p1) / denom^3
  p.rhophi <- s0 * s1 * (s0 * (1 - 2 * p1) - s1 * (1 - 2 * p0)) / denom^3
  p.phiphi <- s0 * s1 * (s0^2 * (1 - 2 * p1) + s1^2 * (1 - 2 * p0)) / denom^3

  ## p_x=p0+x*rho; include BOTH terms in the second-order chain rule.
  list(p = p, d.a = (p.rho + x) * rho.first, d.b = p.phi,
       d.aa = p.rhorho * rho.first^2 + (p.rho + x) * rho.second,
       d.ab = p.rhophi * rho.first, d.bb = p.phiphi,
       interior = all(is.finite(probs) & probs > 1e-8 & probs < 1 - 1e-8))
}

max.likelihood.firth.rd <- function(param, y, x, va, vb, alpha.start, beta.start,
                                    weight, max.step, thres, pa, pb) {
  max.likelihood.firth(y, x, va, vb, alpha.start, beta.start, weight,
                       max.step, thres, pa, pb, firth_rd_derivatives)
}
