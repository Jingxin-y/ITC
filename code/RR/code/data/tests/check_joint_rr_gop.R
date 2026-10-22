## 在 anchored_v1 运行：Rscript data/tests/check_joint_rr_gop.R
source("load_anchored.R", encoding = "UTF-8")
source_anchored_data()

make_design <- function(rho = .5, lambda = .3, ess_ratio = .5, q = NULL) {
  joint_rr_gop_design(alpha_AC = .4, alpha_BC = .2, lambda = lambda,
                      rho = rho, beta0 = -3, beta_E = .2, beta_P = .4,
                      ess_ratio = ess_ratio, q = q)
}
design <- make_design()
tab <- design$risk_table

## 独立枚举 64 个二元组合，核对目标分布的分数矩与边际真值。
z <- as.matrix(expand.grid(rep(list(0:1), 6)))
scores <- rr_gop_scores(z)
row_prob <- function(prob) {
  apply(sweep(z, 2, prob, function(x, p) p^x * (1-p)^(1-x)), 1, prod)
}
mass <- row_prob(design$e)
stopifnot(abs(sum(mass) - 1) < 1e-12)
for (score in c("S_E", "S_P")) {
  stopifnot(abs(sum(mass * scores[[score]]) - .9) < 1e-12,
            abs(sum(mass * (scores[[score]] - .9)^2) - .63) < 1e-12)
}
stopifnot(abs(sum(mass * scores$h_E * scores$h_P)) < 1e-12)
index <- 1 + scores$S_E + 4*scores$S_P
for (population in c("AC", "BC")) {
  prob <- if (population == "AC") design$q else design$e
  risk <- colSums(tab[index, c("p_A", "p_B", "p_C")] * row_prob(prob))
  stopifnot(max(abs(risk - design$marginal_truth[population, names(risk)])) < 1e-12)
}
cat("PASS: fixed score moments and exact marginal risks from 64 covariate patterns\n")

## 直接代回用户给定的三条方程与 c 的可行区间，涵盖三个 rho 和罕见/常见风险。
for (rho in c(0, .5, 1)) for (beta0 in c(-12, -3, 3)) {
  t <- rr_gop_risk_table(-.2, .4, .5, rho, beta0, -.3, .6)
  p <- as.matrix(t[c("p_A", "p_B", "p_C")])
  RA <- exp(-.2 + .5*t$h_E)
  RB <- exp(.4 + rho*.5*t$h_E)
  G <- exp(beta0 - .3*t$h_E + .6*t$h_P)
  c0 <- t$p_C
  lhs <- RA*RB*c0^3 / ((1-RA*c0)*(1-RB*c0)*(1-c0))
  stopifnot(nrow(t) == 16L, all(p > 0 & p < 1),
            all(c0 < 1/pmax(1, RA, RB)),
            max(abs(log(t$p_A/t$p_C) - log(RA))) < 1e-10,
            max(abs(log(t$p_B/t$p_C) - log(RB))) < 1e-10,
            max(abs(log(lhs) - log(G))) < 1e-9)
}
stopifnot(max(abs(make_design(rho=0)$marginal_truth[, "log_RR_BC"] - .2)) < 1e-10,
          max(abs(make_design(rho=1)$marginal_truth[, "log_RR_AB"] - .2)) < 1e-10,
          max(abs(make_design(lambda=0)$marginal_truth[, "log_RR_AC"] - .4)) < 1e-10)
cat("PASS: joint RR/GOP equations, valid risks, rho=0/0.5/1 and lambda=0\n")

## 从 64 格的重要性权重直接计算 ESS，核对反解时 AC/BC 的方向。
for (r in c(1, .8, .5, .25, .1)) {
  q <- ipd_prob_for_ess_theory(r)
  source.mass <- row_prob(q)
  w <- mass / source.mass
  direct.ess <- sum(source.mass*w)^2 / sum(source.mass*w^2)
  stopifnot(abs(direct.ess-r) < 1e-10)
}
q.direct <- c(.2,.4,.6,.8,.5,.3)
direct <- make_design(q=q.direct)
stopifnot(identical(unname(direct$q), q.direct),
          identical(direct$risk_table, design$risk_table),
          identical(make_design(ess_ratio=1)$q, design$e))
q.custom <- ipd_prob_for_ess_theory(.5, direction=c(1,0,-1,0,1,0))
stopifnot(abs(rel_ess_binary(q.custom, design$e)-.5) < 1e-10)
cat("PASS: fixed BC distribution, direct q and theoretical ESS calibration\n")

## 固定 seed 的抽样检验：独立协变量、两组人数固定相等的随机分配、Bernoulli 结局。
set.seed(518)
sim <- generate_anchored_data(50000, 60000, design)
for (name in c("data.IPD", "data.AD.or")) {
  d <- sim[[name]]
  covariates <- covariate_matrix(d$data)
  stopifnot(ncol(covariates) == 6L,
            max(abs(colMeans(covariates) - d$design$covariate_prob)) < .01,
            max(abs(cov(covariates)[upper.tri(cov(covariates))])) < .005,
            sum(d$data$x == 0) == nrow(d$data) / 2,
            sum(d$data$x == 1) == nrow(d$data) / 2,
            all(d$count[c("Na0", "Na1")] == nrow(d$data) / 2),
            max(abs(cov(d$data$x, covariates))) < .005,
            all(d$data$h_E == (d$data$S_E-.9)/sqrt(.63)),
            all(d$data$h_P == (d$data$S_P-.9)/sqrt(.63)))
  idx <- d$truth$risk_table_row
  active <- if (name == "data.IPD") tab$p_A else tab$p_B
  expected <- ifelse(d$data$x == 1, active[idx], tab$p_C[idx])
  stopifnot(identical(d$truth$p0, tab$p_C[idx]),
            identical(d$truth$p1, active[idx]),
            identical(d$truth$observed_risk, expected))
  for (arm in 0:1) {
    take <- d$data$x == arm
    residual <- sum(d$data$y[take] - expected[take])
    stopifnot(abs(residual) < 5*sqrt(sum(expected[take]*(1-expected[take]))))
  }
}
stopifnot(identical(sim$data.AD$anchored.mean.x, colMeans(covariate_matrix(sim$data.AD.or$data))))
set.seed(22)
small <- generate_anchored_data(2, 2, design)
stopifnot(ncol(covariate_matrix(small$data.IPD$data)) == 6L,
          all(is.finite(small$data.AD$anchored.mean.x)),
          all(small$data.IPD$count[c("Na0", "Na1")] == 1),
          all(small$data.AD.or$count[c("Na0", "Na1")] == 1))
for (group in c("IPD", "AD")) for (n in c(1, 3)) {
  stopifnot(inherits(try(data.generation(group, n, design), silent=TRUE), "try-error"))
}
set.seed(23)
a <- generate_anchored_data(100, 150, design)
set.seed(23)
b <- generate_anchored_data(100, 150, design)
stopifnot(identical(a,b),
          all(a$data.IPD$count[c("Na0", "Na1")] == 50),
          all(a$data.AD.or$count[c("Na0", "Na1")] == 75))
cat("PASS: common comparator table, exactly balanced random allocation, AD aggregation and seed reproducibility\n")
