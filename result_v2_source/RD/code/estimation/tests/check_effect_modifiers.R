## 从 anchored_v1 运行：Rscript estimation/tests/check_effect_modifiers.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

## RR/RD 导数和修正项的独立检验见 check_firth_rr_rd.R。

## 非随机的分层计数，用两个 alpha 参数拟合有修饰效应的 RR。
dat <- expand.grid(id = 1:100, x = 0:1, z1 = 0:1)
rates <- c(.2, .3, .3, .6)
dat$y <- as.integer(dat$id <= 100 * rates[1 + dat$x + 2 * dat$z1])
ipd <- list(data = dat)
target <- c(z1 = .4)
wf <- get.weight(covariate_matrix(dat), target)
fit.args <- list(group = "IPD", data = ipd, weight = wf$weight,
  gamma.ipd = wf$ga, target_mean = target, method = "firth", effect_modifiers = "z1")
fit <- suppressWarnings(do.call(get.estimate, fit.args))
stopifnot(ncol(fit$va) == 2L, length(fit$point.est) == 4L,
          identical(dim(fit$cov), c(4L, 4L)), all(is.finite(fit$coefficients)))
stopifnot(identical(unname(fit$va), unname(cbind(1, dat$z1))),
          identical(unname(fit$vb), unname(cbind(1, dat$z1))))
stopifnot(isTRUE(all.equal(as.vector(fit$param.est),
  exp(fit$point.est[1] + dat$z1 * fit$point.est[2]))))

fit.args$method <- "brm"
fit.brm <- suppressWarnings(do.call(get.estimate, fit.args))
truth <- c(log(.3 / .2), log(.6 / .3) - log(.3 / .2))
stopifnot(max(abs(fit.brm$point.est[1:2] - truth)) < .01)
cat("PASS: BRM/Firth fit with modifier, full covariance, and conditional RR predictions\n")

## 没有指定修饰变量时，保留原来的截距-only alpha 设计。
fit.args$effect_modifiers <- NULL
fit.base <- suppressWarnings(do.call(get.estimate, fit.args))
fit.args$effect_modifiers <- character()
fit.empty <- suppressWarnings(do.call(get.estimate, fit.args))
stopifnot(ncol(fit.base$va) == 1L, identical(fit.base, fit.empty))
stopifnot(!grepl("x:", paste(deparse(weighted_glm_formula(dat)), collapse = ""), fixed = TRUE))
cat("PASS: default model and current GLM formulas still exclude effect modifiers\n")

edge <- dat
edge$y[edge$x == 0] <- 0
stopifnot(inherits(try(get.brm.ad(fit, edge, wf$weight), silent = TRUE), "try-error"))
ad <- get.aggre(ipd)
ad$anchored.mean.x <- target
bc.fits <- get.brm.bc(ipd, ad, wf$weight, wf$ga, fit.brm, fit.brm, run_exact=FALSE)
stopifnot(ncol(bc.fits$firth$va) == 2L, length(bc.fits$firth.bc$low) == 2L,
          all(is.na(bc.fits$bayes.bc$low)))
cat("PASS: BRM/Firth BC wrapper preserves modifiers; scalar Bayes correction remains unsupported\n")

## 只检验 pa>1 的 BC 接口：保留上游 bootstrap 次数，缩短优化和搜索；不检验覆盖率。
set.seed(9)
tiny <- dat[seq(1, nrow(dat), by = 25), ]
ci <- exact("RR", tiny$y, tiny$x, cbind(1, tiny$z1), matrix(1, nrow(tiny), 1),
  weight = rep(1, nrow(tiny)), pars = c(.2, .3, -1), se = c(.5, .5, .5),
  pa = 2, pb = 1, max.step = 1, optim.maxit = 1, bisection.max.step = 1,
  grid.mult = 25, thres.dicho = .5, time.limit.sec = 300)
stopifnot(length(ci$low) == 2L, length(ci$up) == 2L, length(ci$p) == 2L,
          all(is.finite(unlist(ci[c("low", "up", "p")]))))
cat("PASS: BC handles both alpha components (small smoke test only)\n")

## 最新 DGM 无需截距数据列；h_E 是 target model 中的一列，须完整传给 BC。
design <- joint_rr_gop_design(.4, .2, .3, .5, -3, .2, .4, ess_ratio=.5)
set.seed(81)
sim <- generate_anchored_data(400, 400, design)
stopifnot(!"v.1" %in% names(sim$data.IPD$data),
          !"v.1" %in% names(sim$data.AD.or$data))
score.fit <- get.estimate("AD", sim$data.IPD, rep(1, 400),
                          method="firth", effect_modifiers="h_E")
stopifnot(identical(unname(score.fit$va), unname(cbind(1, sim$data.IPD$data$h_E))),
          ncol(score.fit$vb) == 4L, length(score.fit$point.est) == 6L,
          all(is.finite(score.fit$coefficients)))
## 上面的 tiny 检查运行真实上游 BC；这里专门检查高维接口原样传递完整矩阵，
## 用记录调用的替身避免再次运行数千次高维优化。替身区间不用于推断。
exact.real <- exact
score.call <- NULL
exact <- function(...) {
  score.call <<- list(...)
  list(low=rep(-1,score.call$pa),up=rep(1,score.call$pa),p=rep(.5,score.call$pa))
}
score.ci <- tryCatch(
  exact_safe("RR", sim$data.IPD$data$y, sim$data.IPD$data$x,
    score.fit$va, score.fit$vb, rep(1, 400), score.fit$point.est, score.fit$se.est,
    ncol(score.fit$va), ncol(score.fit$vb), label="h_E Firth BC",
    bc_control=list(parallel.bootstrap=TRUE,ncores.bootstrap=2L,time.limit.sec=14400)),
  finally={exact <- exact.real}
)
stopifnot(identical(score.call$va,score.fit$va),
          identical(score.call$vb,score.fit$vb),
          identical(score.call$pars,score.fit$point.est),
          identical(score.call$se,score.fit$se.est),
          identical(score.call$y,sim$data.IPD$data$y),
          identical(score.call$x,sim$data.IPD$data$x),
          identical(score.call$weight,rep(1,400)),score.call$pa==2L,score.call$pb==4L,
          isTRUE(score.call$parallel.bootstrap),score.call$ncores.bootstrap==2L,
          score.call$time.limit.sec==14400,
          identical(score.ci,list(low=rep(-1,2),up=rep(1,2),p=rep(.5,2))))
cat("PASS: high-dimensional fitted matrices, parameters and upstream controls forwarded unchanged (call spy)\n")
