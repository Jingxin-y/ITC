## Coefficient-BC compatibility; marginal integration is checked in check_marginal_pipeline.R.
## 从 anchored_v1 运行：Rscript estimation/tests/check_marginal_op.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

## 独立 oracle：16 格的已知概率作为 y=0/1 两行的频数权重。
## 不从拟合模型生成 oracle；先对两臂风险积分，再检验 RR 与 OP。
design <- joint_rr_gop_design(.4, .2, .3, .5, -3, .2, .4, ess_ratio = .5)
tab <- design$risk_table
dat <- merge(tab, expand.grid(x = 0:1, y = 0:1), by = NULL)
risk <- ifelse(dat$x == 1, dat$p_A, dat$p_C)
make_weight <- function(prob_p) {
  10000 * dbinom(dat$S_E, 3, .3) * dbinom(dat$S_P, 3, prob_p) *
    ifelse(dat$y == 1, risk, 1 - risk)
}
w <- make_weight(.3)
args <- list(group = "AD", data = list(data = dat), weight = w,
             effect_modifiers = "h_E", non_effect_modifiers = "h_P",
             model_covariates = c("h_E", "h_P"))
fit <- do.call(get.estimate, args)
stopifnot(ncol(fit$va) == 2L, ncol(fit$vb) == 4L,
          length(fit$point.est) == 6L, all(is.finite(fit$cov)),
          max(abs(fit$point.est[1:2] - c(.4, .3))) < 1e-3,
          identical(fit$op_spec$non_effect_modifiers, "h_P"),
          identical(fit$op_spec$conditioning_covariates, "h_E"))
expected <- sapply(0:3, function(s) {
  rows <- tab$S_E == s
  colSums(tab[rows, c("p_C", "p_A")] * dbinom(tab$S_P[rows], 3, .3))
})
first <- match(0:3, dat$S_E)
pred <- getProbRR(fit$va[first, , drop=FALSE] %*% fit$point.est[1:2],
                 fit$vb[first, , drop=FALSE] %*% fit$point.est[3:6])
stopifnot(max(abs(pred - t(expected))) < 3e-4)
cat("PASS: four-cell reduced BRM recovers integrated risks and original RR parameters\n")

## 非 EM 分布改变：RR 不变，reduced OP 改变；未用结果/真值选择 EM。
args$weight <- make_weight(.65)
changed <- do.call(get.estimate, args)
stopifnot(max(abs(changed$point.est[1:2] - c(.4, .3))) < 1e-3,
          max(abs(changed$point.est[3:6] - fit$point.est[3:6])) > .01)
args$weight <- w
swapped <- brm_op_design(dat, "h_P", "h_E", c("h_E", "h_P"))
stopifnot(identical(swapped$spec$conditioning_covariates, "h_P"),
          identical(as.vector(swapped$va[, 2]), dat$h_P))
cat("PASS: non-EM mixing changes OP, and caller-supplied EM/non-EM roles are respected\n")

## 没有 EM：普通 MLE 必须等于两臂加权事件率之比。
args$effect_modifiers <- character(0)
args$non_effect_modifiers <- c("h_E", "h_P")
marginal <- do.call(get.estimate, args)
p <- vapply(0:1, function(a) weighted.mean(dat$y[dat$x == a], w[dat$x == a]), numeric(1))
stopifnot(ncol(marginal$va) == 1L, ncol(marginal$vb) == 1L,
          abs(marginal$point.est[1] - log(p[2]/p[1])) < 1e-3)
cat("PASS: no EM gives fully marginal weighted risk ratio\n")

## 多个二元 EM 用联合格而非主效应；数据行顺序不改变参考格。
binary <- expand.grid(z1=0:1, z2=0:1, z3=0:1, x=0:1, id=1:10)
binary$y <- as.integer(binary$id <= 2 + binary$x + binary$z1 + binary$z3)
joint <- brm_op_design(binary, c("z1", "z2"), "z3")
stopifnot(ncol(joint$vb) == 4L, ncol(joint$va) == 3L,
          qr(joint$vb)$rank == 4L)
rev_rows <- rev(seq_len(nrow(binary)))
reordered <- brm_op_design(binary[rev_rows, ], c("z1", "z2"), "z3")
stopifnot(identical(joint$vb[rev_rows, ], reordered$vb))
## 未出现的 EM 组合不添加全零 nuisance 列。
partial <- brm_op_design(binary[!(binary$z1 == 1 & binary$z2 == 1), ],
                         c("z1", "z2"), "z3")
stopifnot(ncol(partial$vb) == 3L)
bad <- function(expr) stopifnot(inherits(try(force(expr), silent=TRUE), "try-error"))
bad(brm_op_design(binary, "z1", "z1"))
bad(brm_op_design(binary, "missing", "z3"))
bad(brm_op_design(binary, c("z1", "z1"), "z3"))
bad(brm_op_design(binary, "x", "z3"))
continuous <- binary
continuous$em <- seq_len(nrow(binary))
bad(get.estimate("AD", list(data=continuous), rep(1,nrow(continuous)), effect_modifiers="em"))
cat("PASS: multiple EM, observed-cell coding, stable reference and invalid-input rejection\n")

## 真实比较流程：所有 BRM/Firth 及 BC 沿用四格设计，其他方法仍用自己的主效应。
ipd <- list(data=binary)
ad <- get.aggre(ipd)
ad$anchored.mean.x <- c(z1=.4,z2=.6,z3=.3)
recorded <- list()
real_exact <- exact
exact <- function(...) {
  a <- list(...)
  recorded[[length(recorded)+1L]] <<- a
  list(low=rep(-1,a$pa),up=rep(1,a$pa),p=rep(.5,a$pa))
}
comparison <- tryCatch(suppressWarnings(get.compare.anchored(bc_target="coefficient",score_standardization="fixed",ipd,ad,
  weight_covariates="z3", model_covariates=c("z1","z2","z3"),
  effect_modifiers=c("z1","z2"), non_effect_modifiers="z3")),
  finally={exact <- real_exact})
for (name in c("brm", "brm-BC", "brm-FC-BC")) {
  f <- comparison$fits[[name]]
  stopifnot(ncol(f$va)==3L,ncol(f$vb)==4L,all(is.finite(f$point.est)),
            identical(f$op_spec$non_effect_modifiers,"z3"))
}
stopifnot(length(recorded)==2L, recorded[[1]]$pa==3L, recorded[[1]]$pb==4L,
          identical(recorded[[1]]$vb,comparison$fits$brm$vb),
          identical(recorded[[2]]$vb,comparison$fits[["brm-FC-BC"]]$vb),
          identical(comparison$variables$op_model,"marginal"),
          "z3" %in% colnames(model.matrix(comparison$fits$LP)))
## 单独调用 BC wrapper 时，必须从 fit 元数据恢复 OP，而不能把 dummy 名称当输入列。
bc <- suppressWarnings(get.brm.bc(bc_target="coefficient",ipd,ad,comparison$weights$weight,comparison$weights$ga,
                                  comparison$fits$brm,NULL,run_exact=FALSE))
stopifnot(identical(bc$firth$vb,comparison$fits$brm$vb))
cat("PASS: comparison, Firth, joint covariance and BC forwarding (interval call spy)\n")
