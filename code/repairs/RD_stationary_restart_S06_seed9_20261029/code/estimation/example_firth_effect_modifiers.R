## 在 anchored_v1 目录运行。独立示例：两个 effect modifier，不改变原批量模拟。
source("load_anchored.R", encoding = "UTF-8")
source_anchored()
set.seed(1)

## 1. 平衡设计；log RR(z) = 0.2 + 0.5*z1 - 0.3*z2。
dat <- expand.grid(id = 1:100, x = 0:1, z1 = 0:1, z2 = 0:1)
va <- cbind(1, dat$z1, dat$z2)
vb <- va
alpha.true <- c(0.2, 0.5, -0.3)
beta.true <- c(-1, 0.4, -0.2)
prob <- getProbRR(va %*% alpha.true, vb %*% beta.true)
dat$y <- rbinom(nrow(dat), 1, prob[cbind(seq_len(nrow(dat)), dat$x + 1)])
data.IPD <- list(data = dat)

## 2. MAIC 权重。本例把目标均值视为已知，方便单独查看 effect modifier。
target.mean <- c(z1 = 0.4, z2 = 0.6)
weight.fit <- get.weight(covariate_matrix(data.IPD$data), target.mean)

## 3. 同一设计分别拟合 BRM 与 Firth；返回完整 alpha、beta 和协方差。
fit.brm <- get.estimate(
  "IPD", data.IPD, weight.fit$weight, weight.fit$ga,
  target_mean = target.mean, method = "brm",
  effect_modifiers = c("z1", "z2")
)
fit.firth <- get.estimate(
  "IPD", data.IPD, weight.fit$weight, weight.fit$ga,
  target_mean = target.mean, method = "firth",
  effect_modifiers = c("z1", "z2")
)

## alpha 1 是 z1=z2=0 时的 log RR；alpha 2、3 为对应的修饰效应。
pa <- ncol(fit.firth$va)
print(cbind(truth = alpha.true,
            brm = fit.brm$point.est[seq_len(pa)],
            firth = fit.firth$point.est[seq_len(pa)]))
print(fit.firth$coefficients[seq_len(pa), , drop = FALSE])
print(c(brm = fit.brm$converged, firth = fit.firth$converged))

## 4. 可选：对每个 alpha 分量计算 BC 区间，而非目标人群的边际 log RR。
## 此处不使用原来的标量 Beta-Binomial 边界修正。
run_bc <- FALSE
if (run_bc) {
  ci.firth <- exact_safe(
    "RR", dat$y, dat$x, fit.firth$va, fit.firth$vb, weight.fit$weight,
    fit.firth$point.est, fit.firth$se.est,
    pa = pa, pb = ncol(fit.firth$vb), label = "Firth with effect modifiers",
    max.step = min(pa * 20, 1000)
  )
  print(cbind(low = ci.firth$low, up = ci.firth$up, p = ci.firth$p))
}
