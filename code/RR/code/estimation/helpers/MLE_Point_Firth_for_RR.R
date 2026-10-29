## RR: theta=log(RR), phi=log(odds product).
## Derivative source: brmplus commit 3b628a6d32dc7fe0f31fbd3ab181aecd3866f982,
## compare/MLE_Point_Firth_for_RR.R; see BRMPLUS_LICENSE.txt.
## Scalar derivatives reuse the original brmplus RR blocks; va/vb enter only
## in estimation/helpers/firth_augmentation.R. Supports effect modifiers (pa > 1).
firth_rr_derivatives <- function(alpha, beta, x, va, vb) {
  probs <- getProbRR(as.vector(va %*% alpha), as.vector(vb %*% beta))
  p0 <- pmin(pmax(probs[, 1], 1e-8), 1 - 1e-8)
  p1 <- pmin(pmax(probs[, 2], 1e-8), 1 - 1e-8)
  p <- ifelse(x == 1, p1, p0)
  q0 <- 1 - p0
  q1 <- 1 - p1
  denom <- q0 + q1

  ## First and second derivatives of log(p0).
  g.a <- -q0 / denom
  g.b <- q0 * q1 / denom
  g.aa <- ((p0 - p1) * g.a - q0 * p1) / denom^2
  g.ab <- q0 * q1 * (p0 - p1) / denom^3
  g.bb <- -(p0 * q1^2 + p1 * q0^2) * g.b / denom^2

  ## Convert to derivatives of p_x = exp(log(p0) + x*theta).
  list(p = p, d.a = p * (g.a + x), d.b = p * g.b,
       d.aa = p * ((g.a + x)^2 + g.aa),
       d.ab = p * ((g.a + x) * g.b + g.ab),
       d.bb = p * (g.b^2 + g.bb),
       interior = all(is.finite(probs) & probs > 1e-8 & probs < 1 - 1e-8))
}

max.likelihood.firth.rr <- function(param, y, x, va, vb, alpha.start, beta.start,
                                    weight, max.step, thres, pa, pb) {
  max.likelihood.firth(y, x, va, vb, alpha.start, beta.start, weight,
                       max.step, thres, pa, pb, firth_rr_derivatives)
}
