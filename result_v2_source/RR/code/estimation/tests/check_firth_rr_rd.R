## 从 anchored_v1 运行：Rscript estimation/tests/check_firth_rr_rd.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

## 独立数值微分：只调用概率映射，不使用被检验的解析导数。
numerical_derivatives <- function(prob, theta, step = 1e-4) {
  p <- length(theta)
  base <- prob(theta)
  gradient <- numeric(p)
  hessian <- matrix(0, p, p)
  for (j in seq_len(p)) {
    ej <- numeric(p); ej[j] <- step
    gradient[j] <- (prob(theta + ej) - prob(theta - ej)) / (2 * step)
    hessian[j, j] <- (prob(theta + ej) - 2 * base + prob(theta - ej)) / step^2
    if (j > 1) for (k in seq_len(j - 1)) {
      ek <- numeric(p); ek[k] <- step
      hessian[j, k] <- hessian[k, j] <-
        (prob(theta + ej + ek) - prob(theta + ej - ek) -
           prob(theta - ej + ek) + prob(theta - ej - ek)) / (4 * step^2)
    }
  }
  list(p = base, gradient = gradient, hessian = hessian)
}

set.seed(19)
n <- 24
x <- rep(0:1, length.out = n)
z <- matrix(runif(n * 3, -.6, .6), n, 3)
vb <- cbind(1, z)
beta <- c(-1, .3, -.2, .1)
weight <- rep(c(.7, 1.3, 2), length.out = n)
y <- rbinom(n, 1, .4)

for (model in c("RR", "RD")) for (pa in 1:3) {
  va <- cbind(1, z[, seq_len(pa - 1), drop = FALSE])
  alpha <- c(.2, rep(.1, pa - 1))
  theta <- c(alpha, beta)
  deriv <- if (model == "RR") firth_rr_derivatives else firth_rd_derivatives
  get_prob <- if (model == "RR") getProbRR else getProbRD
  terms <- deriv(alpha, beta, x, va, vb)
  got <- firth_score_adjustment(y, va, vb, weight, terms)
  a <- seq_len(pa)
  b <- pa + seq_len(ncol(vb))
  Iinv <- solve(got$information)
  V <- got$covariance
  oracle <- numeric(length(theta))
  for (i in seq_len(n)) {
    probability <- function(par) {
      as.numeric(get_prob(sum(va[i, ] * par[a]), sum(vb[i, ] * par[b]))[1, x[i] + 1])
    }
    numeric <- numerical_derivatives(probability, theta)
    D <- c(va[i, ] * terms$d.a[i], vb[i, ] * terms$d.b[i])
    H <- rbind(
      cbind(tcrossprod(va[i, ]) * terms$d.aa[i],
            outer(va[i, ], vb[i, ]) * terms$d.ab[i]),
      cbind(outer(vb[i, ], va[i, ]) * terms$d.ab[i],
            tcrossprod(vb[i, ]) * terms$d.bb[i]))
    stopifnot(max(abs(numeric$gradient - D)) < 2e-6,
              max(abs(numeric$hessian - H)) < 2e-5)
    oracle <- oracle + .5 * weight[i] * numeric$gradient /
      (numeric$p * (1 - numeric$p)) * sum(diag(V %*% numeric$hessian))
    variance <- numeric$p * (1 - numeric$p)
    C <- numeric$hessian / variance -
      (1 - 2 * numeric$p) * tcrossprod(numeric$gradient) / variance^2
    oracle <- oracle + as.vector(weight[i] * C %*%
      (V - weight[i] * Iinv) %*% numeric$gradient)
  }
  stopifnot(max(abs(got$adjustment - oracle)) < 1e-4)
}
cat("PASS: RR/RD gradients, Hessians and fixed-weight adjustment for pa=1/2/3\n")

## RD: 在 rho=0、phi=0 处检查解析导数，防止遗漏 tanh 的链式法则。
for (theta in list(c(0, 0), c(.2, 0), c(0, -.7))) for (arm in 0:1) {
  terms <- firth_rd_derivatives(theta[1], theta[2], arm, matrix(1,1,1),matrix(1,1,1))
  numeric <- numerical_derivatives(function(par) getProbRD(par[1],par[2])[1,arm+1],
                                  theta, step=1e-3)
  H <- matrix(c(terms$d.aa,terms$d.ab,terms$d.ab,terms$d.bb),2,2)
  stopifnot(max(abs(numeric$gradient-c(terms$d.a,terms$d.b)))<1e-5,
            max(abs(numeric$hessian-H))<1e-4)
}
cat("PASS: RD derivatives at zero RD / zero log odds-product\n")

## 固定分层风险，RD 在 z1=0/1 时为 .1/.3，含一个 effect modifier。
dat <- expand.grid(id = 1:100, x = 0:1, z1 = 0:1)
rates <- c(.2, .3, .3, .6)
dat$y <- as.integer(dat$id <= 100 * rates[1 + dat$x + 2 * dat$z1])
ipd <- list(data = dat)
wf <- get.weight(covariate_matrix(dat), c(z1=.4))
args <- list(group="IPD",data=ipd,weight=wf$weight,gamma.ipd=wf$ga,
             target_mean=.4,method="firth",effect_modifiers="z1",param="RD")
fit <- suppressWarnings(do.call(get.estimate,args))
stopifnot(all(is.finite(fit$coefficients)), ncol(fit$va)==2L,
  all(abs(fit$param.est)<1),isTRUE(all.equal(as.vector(fit$param.est),
  tanh(fit$point.est[1]+dat$z1*fit$point.est[2]))))
args$method <- "brm"
mle <- suppressWarnings(do.call(get.estimate,args))
stopifnot(max(abs(as.vector(mle$param.est)-ifelse(dat$z1==0,.1,.3)))<.01)
args$method <- "firth"
args$effect_modifiers <- NULL
constant <- suppressWarnings(do.call(get.estimate,args))
stopifnot(all(is.finite(constant$coefficients)),ncol(constant$va)==1L)
cat("PASS: RD MLE/Firth via get.estimate(), with/without modifiers, RD-scale predictions\n")

## 固定代表性权重整体缩放不改变修正分数的根。
v <- cbind(1, dat$z1)
w <- rep(c(1,2),length.out=nrow(dat))
fc.args <- list(param="RD",y=dat$y,x=dat$x,va=v,vb=v,weight=w,
               alpha.start=c(0,0),beta.start=c(0,0),pa=2,pb=2,max.step=40,thres=1e-6)
weighted <- suppressWarnings(do.call(max.likelihood.firth.rd,fc.args))
fc.args$weight <- 10 * w
scaled <- suppressWarnings(do.call(max.likelihood.firth.rd,fc.args))
stopifnot(isTRUE(weighted$convergence), isTRUE(scaled$convergence),
          max(abs(weighted$par-scaled$par))<1e-8)
cat("PASS: RD fixed analytic weights are invariant to global rescaling\n")

set.seed(29)
tiny <- dat[seq(1,nrow(dat),by=25),]
ci <- exact("RD",tiny$y,tiny$x,cbind(1,tiny$z1),matrix(1,nrow(tiny),1),
  weight=rep(1,nrow(tiny)),pars=c(.1,.15,-1),se=c(.2,.2,.5),
  pa=2,pb=1,max.step=1,optim.maxit=1,bisection.max.step=1,
  grid.mult=25,thres.dicho=.5,time.limit.sec=300)
stopifnot(length(ci$low)==2L,length(ci$up)==2L,length(ci$p)==2L,
          all(is.finite(unlist(ci[c("low","up","p")]))))
cat("PASS: RD BC interface with multiple alpha parameters (small smoke test)\n")
