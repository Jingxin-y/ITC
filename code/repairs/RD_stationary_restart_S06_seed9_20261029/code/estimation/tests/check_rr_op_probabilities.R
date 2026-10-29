## Run from anchored_v1: Rscript estimation/tests/check_rr_op_probabilities.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
check <- function(ok, label) { stopifnot(ok); cat("PASS:", label, "\n") }
expect_error <- function(expr) inherits(try(force(expr), silent=TRUE), "try-error")

## Independent oracle: solve the OP equation in the larger risk's logit,
## rather than reuse any quadratic formula or implementation helper.
oracle <- function(t, s) {
  a <- abs(t)
  equation <- function(u) {
    log.low <- plogis(u, log.p=TRUE) - a
    u + log.low - log(-expm1(log.low)) - s
  }
  u <- uniroot(equation, c(-2000, 2000), tol=1e-11)$root
  log.high <- plogis(u, log.p=TRUE)
  probs <- exp(c(log.high-a, log.high))
  if (t < 0) rev(probs) else probs
}
grid <- expand.grid(t=c(-30,-12.0001,-12,-5,-1e-8,0,1e-8,5,12,12.0001,30),
                    s=c(-60,-30,-13,-12.0001,-12,-1e-8,0,1e-8,12,12.0001,30,60))
p <- getProbRR(grid$t, grid$s)
ref <- t(vapply(seq_len(nrow(grid)), function(i) oracle(grid$t[i],grid$s[i]), numeric(2)))
check(max(abs(p-ref)) < 1e-11, "independent OP root oracle across/beyond legacy +/-12 thresholds")
interior <- apply(p, 1, function(z) all(z > 1e-12 & z < 1-1e-10))
check(max(abs(log(p[interior,2])-log(p[interior,1])-grid$t[interior])) < 1e-12,
      "defining log RR identity")
check(max(abs(qlogis(p[interior,1])+qlogis(p[interior,2])-grid$s[interior])) < 2e-6,
      "defining log OP identity in representable interior")

## Closed forms independent of the general implementation.
values <- c(-1300,-60,-13,-12,-1e-10,0,1e-10,12,30,60)
equal <- getProbRR(rep(0,length(values)),values)
check(max(abs(equal[,1]-plogis(values/2))) < 1e-14 &&
      identical(equal[,1],equal[,2]), "RR=1 gives equal risks expit(logOP/2), including logOP=-13")
tiny <- getProbRR(0,-1300)[1,1]
check(tiny > 0 && abs(log(tiny)+650) < 1e-12,
      "representable tiny risks survive intermediate exponential overflow")
zero.op <- getProbRR(values,rep(0,length(values)))
check(max(abs(zero.op[,1]-plogis(-values))) < 1e-14 &&
      max(abs(zero.op[,2]-plogis(values))) < 1e-14, "OP=1 closed form")
swapped <- getProbRR(-grid$t,grid$s)
check(max(abs(swapped-p[,2:1])) < 1e-14, "swapping treatment arms preserves OP and inverts RR")

## A small parameter perturbation must not jump to an artificial 0/1.
for (boundary in c(-12,12)) {
  op <- getProbRR(rep(log(2),3), boundary+c(-1e-6,0,1e-6))
  rr <- getProbRR(boundary+c(-1e-6,0,1e-6),rep(-2,3))
  check(all(op > 0 & op < 1) && all(rr > 0 & rr < 1) &&
        max(abs(op[1,]-op[3,])) < 1e-6 && max(abs(rr[1,]-rr[3,])) < 1e-6,
        paste("continuity across log RR/log OP =",boundary))
}
extreme <- getProbRR(c(-1000,0,1000,0,-1000,1000),c(-1000,-1300,-1000,1300,1000,1000))
check(all(is.finite(extreme)) && all(extreme >= 0 & extreme <= 1),
      "extreme finite predictors return finite representable probabilities")

check(identical(getProbScalarRR(c(log(2),-13)), as.vector(getProbRR(log(2),-13))) &&
      identical(getProbRR(c(log(2),-13)),getProbRR(matrix(c(log(2),-13),1,2))) &&
      identical(getProbRR(0,c(-13,0)),getProbRR(c(0,0),c(-13,0))) &&
      identical(dim(getProbRR(numeric(),numeric())),c(0L,2L)),
      "scalar, packed-pair, matrix, scalar recycling and empty interfaces")
check(expect_error(getProbRR(c(1,2),c(1,2,3))) && expect_error(getProbRR(Inf,0)) &&
      expect_error(getProbScalarRR(0)), "invalid inputs fail explicitly")

## Analytic FC derivatives must remain valid past the former threshold.
for (theta in list(c(0,-13),c(log(2),-13),c(12.01,-2),c(-12.01,-2))) for (arm in 0:1) {
  f <- function(z) getProbScalarRR(z)[arm+1L]
  h <- 1e-4; e1 <- c(h,0); e2 <- c(0,h); base <- f(theta)
  grad <- c((f(theta+e1)-f(theta-e1))/(2*h),(f(theta+e2)-f(theta-e2))/(2*h))
  cross <- (f(theta+e1+e2)-f(theta+e1-e2)-f(theta-e1+e2)+f(theta-e1-e2))/(4*h^2)
  H <- matrix(c((f(theta+e1)-2*base+f(theta-e1))/h^2,cross,cross,
                (f(theta+e2)-2*base+f(theta-e2))/h^2),2,2)
  z <- firth_rr_derivatives(theta[1],theta[2],arm,matrix(1,1,1),matrix(1,1,1))
  check(isTRUE(z$interior) && max(abs(grad-c(z$d.a,z$d.b))) < 1e-7 &&
        max(abs(H-matrix(c(z$d.aa,z$d.ab,z$d.ab,z$d.bb),2,2))) < 1e-7,
        paste("FC derivative oracle at",paste(theta,collapse=","),"arm",arm))
}
check(!firth_rr_derivatives(0,-1300,0,matrix(1,1,1),matrix(1,1,1))$interior,
      "FC retains the original interior-probability safeguard")
