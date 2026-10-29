## Mean-bias reduction for independent Bernoulli outcomes with fixed analytic
## (including MAIC) weights. See ../FC_FIXED_WEIGHTS.md for the derivation.
## I = sum w DD'/v, B = sum w^2 DD'/v, V = I^{-1} B I^{-T}.
## A = 1/2 sum w D/v trace(V H) + sum w C (V - w I^{-1}) D,
## where C = H/v - (1-2p) DD'/v^2. Unit weights recover ordinary Firth.
firth_score_adjustment <- function(y, va, vb, weight, terms) {
  pa <- ncol(va)
  a <- seq_len(pa)
  b <- pa + seq_len(ncol(vb))
  variance <- terms$p * (1 - terms$p)
  D <- cbind(va * terms$d.a, vb * terms$d.b)
  information <- crossprod(D, D * (weight / variance))
  information.inv <- solve(information)
  score.variance <- crossprod(D, D * (weight^2 / variance))
  covariance <- information.inv %*% score.variance %*% t(information.inv)
  covariance <- (covariance + t(covariance)) / 2

  ## trace(V H_i) from three blocks; no n*p*p array or cumulant tensors.
  h.aa <- rowSums((va %*% covariance[a, a, drop = FALSE]) * va)
  h.ab <- rowSums((va %*% covariance[a, b, drop = FALSE]) * vb)
  h.bb <- rowSums((vb %*% covariance[b, b, drop = FALSE]) * vb)
  curvature <- terms$d.aa * h.aa + 2 * terms$d.ab * h.ab + terms$d.bb * h.bb

  ## Each row of R is ((V - w_i I^{-1}) D_i)'. Contract H_i R_i
  ## using its alpha/alpha, alpha/beta and beta/beta blocks.
  R <- D %*% covariance - weight * (D %*% information.inv)
  r.a <- rowSums(va * R[, a, drop = FALSE])
  r.b <- rowSums(vb * R[, b, drop = FALSE])
  HR <- cbind(va * (terms$d.aa * r.a + terms$d.ab * r.b),
              vb * (terms$d.ab * r.a + terms$d.bb * r.b))
  CR <- HR / variance - D * ((1 - 2 * terms$p) * rowSums(D * R) / variance^2)

  score <- colSums(D * (weight * (y - terms$p) / variance))
  curvature.adjustment <- 0.5 * colSums(D * (weight * curvature / variance))
  weight.adjustment <- colSums(CR * weight)
  adjustment <- curvature.adjustment + weight.adjustment
  list(score = score, adjustment = adjustment,
       adjusted.score = score + adjustment, information = information,
       score.variance = score.variance, covariance = covariance,
       curvature.adjustment = curvature.adjustment,
       weight.adjustment = weight.adjustment)
}
