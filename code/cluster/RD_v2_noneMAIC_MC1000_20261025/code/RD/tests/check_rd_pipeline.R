## Base R/MASS only. Tests write only to RD_TEST_OUTPUT or a temporary directory.
args <- commandArgs(trailingOnly = TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
code_dir <- if (length(args)) args[1] else dirname(rd_dir)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
source_anchored_rd(code_dir, rd_dir)
out <- Sys.getenv("RD_TEST_OUTPUT", tempfile("rd_checks_")); dir.create(out, recursive = TRUE)
check <- function(condition, message) { if (!isTRUE(condition)) stop(message); cat("PASS:", message, "\n") }
fails <- function(expression) inherits(tryCatch(expression, error = identity), "error")
before <- rd_source_hashes(code_dir, rd_dir)
for (f in list.files(rd_dir, "\\.R$", recursive = TRUE, full.names = TRUE)) parse(f, encoding = "UTF-8")
check(TRUE, "all RD R sources parse")

cal <- rd_calibration_table()
adjusted_difference <- cal$difference - cal$intended_risk_shift
check(nrow(cal) == 360 && max(abs(adjusted_difference[cal$population == "BC"])) < 1e-11,
  "60 scenarios preserve BC C/B means and apply the intended A-risk increase")
check(max(abs(adjusted_difference[cal$population == "AC" & cal$event == "common"])) < .066 &&
  max(abs(adjusted_difference[cal$population == "AC" & cal$event == "rare"])) < .0077,
  "AC-RR discrepancies after the explicit A shift retain v1 calibrated bounds")
write.csv(cal, file.path(out, "rd_rr_calibration_comparison.csv"), row.names = FALSE)
truth <- rd_population_truth_table()
check(nrow(truth) == 120 && all(truth$preset_version == 2L) &&
  max(abs(truth$RD_AB[truth$population == "BC" & truth$event == "common"] - .14)) < 1e-12 &&
  max(abs(truth$RD_AB[truth$population == "BC" & truth$event == "rare"] - .05)) < 1e-12 &&
  all(truth$RD_AB[truth$event == "rare"] >= .05 - 1e-12),
  "v2 population truths meet BC A-B targets and rare A-B RD is at least five percentage points")
default_truth <- truth[truth$ess_ratio == .5 & truth$ess_source == "diffuse", ]
check(all(default_truth$RD_AB[default_truth$event == "common"] >= .10 &
  default_truth$RD_AB[default_truth$event == "common"] <= .20),
  "default diffuse/ESS=.5 common A-B RD is within ten to twenty percentage points")
write.csv(truth, file.path(out, "rd_population_truth.csv"), row.names = FALSE)
raw <- expand.grid(rep(list(0:1), 6)); names(raw) <- paste0("z", 1:6)
E <- 1 - raw$z1 + raw$z2 + 1 - raw$z3; P <- raw$z4 + 1 - raw$z5 + raw$z6
for (event in c("common", "rare")) for (src in c("diffuse", "em_driven", "nonem_driven"))
  for (ess in c(.15, .25, .5, .8, 1)) for (level in c("low", "high")) {
    d <- webster_rd_scenario_design(event, ess, src, level); risk <- d$risk_table[1 + E + 4 * P, ]
    check_risk <- as.matrix(risk[c("p_C", "p_A", "p_B")])
    stopifnot(all(check_risk > 0 & check_risk < 1),
      all(check_risk >= if (event == "common") .05 else .005),
      all(check_risk <= if (event == "common") .95 else .995),
      max(abs((risk$p_A - risk$p_C) - (d$parameters$rd_AC + d$parameters$em_rd_AC * E))) < 1e-12,
      max(abs((risk$p_B - risk$p_C) - (d$parameters$rd_BC + d$parameters$em_rd_BC * E))) < 1e-12,
      abs(d$scenario$lambda_BC / d$scenario$lambda_AC - .5) < 1e-12)
  }
check(TRUE, "all 64 covariate cells valid; non-EMs cancel from RD; BC/AC modifier ratio is one half")
for (event in c("common", "rare")) {
  lo <- webster_rd_scenario_design(event, em_level = "low")
  hi <- webster_rd_scenario_design(event, em_level = "high")
  stopifnot(abs(hi$scenario$lambda_AC / lo$scenario$lambda_AC - 4) < 1e-12)
}
check(fails(webster_rd_scenario_design(parameters = list(baseline_risk = -2))) &&
  fails(generate_webster_rd_scenario(101)), "invalid risks and odd sample sizes are rejected")
a <- generate_webster_rd_scenario(500, seed = 19); b <- generate_webster_rd_scenario(500, seed = 19)
check(identical(a, b), "same seed reproduces both randomized trials")

## Saturated effects exactly express additive RD; beta target gradient is zero.
va <- cbind(1, as.numeric(E == 1), as.numeric(E == 2), as.numeric(E == 3))
vb <- va; d <- webster_rd_scenario_design()
theta <- atanh(d$parameters$rd_AC + d$parameters$em_rd_AC * 0:3)
alpha <- c(theta[1], theta[-1] - theta[1]); pars <- c(alpha, -.5, .2, .1, -.2)
target <- marginal_rd_target(pars, va, vb, gradient = TRUE)
numeric_gradient <- vapply(seq_along(pars), function(j) {
  up <- down <- pars; up[j] <- up[j] + 1e-5; down[j] <- down[j] - 1e-5
  (marginal_rd_target(up, va, vb) - marginal_rd_target(down, va, vb)) / 2e-5
}, numeric(1))
check(max(abs(tanh(drop(va %*% alpha)) - (d$parameters$rd_AC + d$parameters$em_rd_AC * E))) < 1e-12 &&
  max(abs(target$gradient - numeric_gradient)) < 1e-8 && all(target$gradient[5:8] == 0),
  "saturated BRM represents true RD; full target gradient agrees with finite differences")
for (th in c(-.8, 0, .8)) for (ph in c(-25, -1, 0, 1, 25)) {
  z <- rd_probabilities(th, ph)
  stopifnot(abs(z$p1 - z$p0 - tanh(th)) < 1e-12, z$p0 >= 0, z$p1 <= 1)
  if (abs(ph) <= 1) stopifnot(max(abs(unlist(z) - getProbScalarRD(th, ph))) < 1e-10)
}
check(TRUE, "stable RD probability map preserves RD at extreme and zero OP")

## Independent polynomial constrained MLE oracle for the MN score.
## Same MN correction as PropCIs::diffscoreci; no shared production optimizer.
mn_oracle <- function(p1, n1, p0, n0, d) {
  if (abs(p1 - p0 - d) < 1e-12) return(0)
  coef <- c(n0 * p0 * d * (1 - d),
    n1 * (p1 - d) + n0 * (p0 * (1 - 2 * d) - d * (1 - d)),
    -n1 * (p1 - d + 1) + n0 * (-p0 - 1 + 2 * d), n1 + n0)
  roots <- polyroot(coef); qs <- Re(roots[abs(Im(roots)) < 1e-8])
  qs <- qs[qs >= max(0, -d) & qs <= min(1, 1 - d)]
  ll <- function(q) n1 * (p1 * log(q + d) + (1 - p1) * log1p(-q - d)) +
    n0 * (p0 * log(q) + (1 - p0) * log1p(-q))
  q <- qs[which.max(vapply(qs, ll, numeric(1)))]
  (p1 - p0 - d)^2 / (((q + d) * (1 - q - d) / n1 + q * (1 - q) / n0) * (n1 + n0) / (n1 + n0 - 1))
}
for (counts in list(c(.2, 100, .1, 100), c(.55, 80, .32, 150), c(.03, 400, .015, 200))) {
  for (dd in c(-.02, 0, .05)) {
    actual <- do.call(rd_mn_score, as.list(c(counts, dd)))
    expected <- do.call(mn_oracle, as.list(c(counts, dd)))
    stopifnot(abs(actual - expected) < 1e-5 * max(1, expected))
  }
  z <- do.call(rd_mn_interval, as.list(counts))
  f <- function(dd) do.call(mn_oracle, as.list(c(counts, dd))) - qchisq(.95, 1)
  stopifnot(abs(f(z$lower)) < 1e-4, abs(f(z$upper)) < 1e-4)
}
check(TRUE, "unweighted MN score/interval agree with independent constrained-polynomial oracle")
for (ps in list(c(0, 0), c(1, 1), c(0, 1), c(1, 0))) {
  z <- rd_mn_interval(ps[1], 50, ps[2], 40)
  stopifnot(is.finite(z$lower), is.finite(z$upper), z$lower <= z$estimate, z$upper >= z$estimate)
}
z <- rd_weighted_arms(a$data.IPD$data$y, a$data.IPD$data$x, rep(1, 500))
stopifnot(z$n0 == 250, z$n1 == 250)
check(fails(rd_mn_interval(.2, .1, .1, .1)), "MN boundary outcomes are supported; invalid effective totals rejected")
weights <- exp(seq(-2, 2, length.out = 500))
z1 <- rd_weighted_arms(a$data.IPD$data$y, a$data.IPD$data$x, weights)
z2 <- rd_weighted_arms(a$data.IPD$data$y, a$data.IPD$data$x, weights * 1e6)
check(max(abs(unlist(z1[1:4]) - unlist(z2[1:4]))) < 1e-10, "effective MN is invariant to weight scaling")

dat <- a$data.IPD$data; gd <- rd_glm_design(dat); lpm <- rd_fit_lpm(dat, weights, gd)
fit <- lm(y ~ x + h_E + h_P + x:h_E, data = dat, weights = weights)
X <- model.matrix(fit); w <- weights / mean(weights); inv <- solve(crossprod(X, X * w))
h <- w * rowSums((X %*% inv) * X)
V <- inv %*% crossprod(X * (w * residuals(fit) / (1 - h))) %*% inv
check(max(abs(lpm$coefficients - coef(fit))) < 1e-10 && max(abs(lpm$covariance - V)) < 1e-9,
  "weighted LPM and HC3 agree with independent lm/matrix calculation")
s <- rd_weighted_sandwich(dat$y, dat$x, weights)
phi <- weights * (dat$x * (dat$y - s$risk_A) / sum(weights * dat$x) -
  (1 - dat$x) * (dat$y - s$risk_C) / sum(weights * (1 - dat$x)))
check(abs(s$se^2 - sum(phi^2)) < 1e-12, "fixed-weight RD sandwich agrees with arm-mean influence variance")

result <- run.anchored.rd(d, 500, seed = 19, run_exact = FALSE)
joint <- run.anchored.rd(d, 500, seed = 19, run_exact = FALSE, model_weighting = "maic")
subset <- run.anchored.rd(d, 500, seed = 19, run_exact = FALSE, weight_covariates = paste0("z", 1:3))
check(identical(result$data$data.IPD$data, joint$data$data.IPD$data) &&
  identical(result$data$data.AD.or$data, subset$data$data.AD.or$data), "analysis weighting does not alter either generated trial")
check(nrow(result$estimates) == 9 && all(result$estimates$status == "available") &&
  all(joint$estimates$status == "available") && ncol(result$fits$brm$va) == 2 &&
  identical(result$variables$brm_effect_model, "linear_atanh") &&
  all(result$estimates$interval_status[c(2, 4)] == "BC_not_run"), "unit and MAIC pipelines return nine RD inference rows with two linear atanh-RD effect parameters")
check(all(joint$estimates$variance_target[c(1, 3, 6, 7, 9)] >= 0), "joint RD variance includes nonnegative AD-mean contribution")
id <- result$fits$IB; ib_ref <- glm(id$design$formula, data = result$data$data.IPD$data,
  family = binomial("identity"), start = c(.4, rep(0, ncol(id$design$X) - 1)), control = glm.control(maxit = 200, epsilon = 1e-9))
check(max(abs(id$coefficients - coef(ib_ref))) < 1e-7 && max(abs(id$covariance - vcov(ib_ref))) < 1e-6,
  "identity-binomial estimates/covariance agree with independent glm")

## Single-cell data gives deterministic successful observed profiling and small BC inversion.
y <- c(rep(1, 20), rep(0, 30), rep(1, 25), rep(0, 25)); x <- rep(0:1, each = 50)
v <- matrix(1, 100, 1); pars <- c(atanh(.1), log(.4 / .6 * .5 / .5))
set.seed(712); rng <- .Random.seed
ci <- exact_marginal_rd_target(y = y, x = x, va = v, vb = v, pars = pars,
  psi.se = .07, psi.range = c(-.1, .3), nsim.coarse = 10L, nsim.confirm = 10L,
  nsim.refine = 10L, nsim.p.value = 20L, thres.dicho = .05,
  bisection.max.step = 8L, seed = 99, time.limit.sec = 60)
check(identical(rng, .Random.seed) && isTRUE(ci$diagnostics$anchor$converged) &&
  ci$diagnostics$counts$bootstrap.draws > 0 && is.finite(ci$p) && isTRUE(ci$interval.available),
  "real small-budget RD BC inversion returns interval/P and preserves caller RNG")
check(abs(ci$estimate - .1) < 1e-12 && ci$effect.estimate == ci$estimate &&
  abs(ci$diagnostics$anchor$constraint.error) < 1e-7,
  "BC outputs are natural RD and observed nuisance fit satisfies target constraint")
saveRDS(ci, file.path(out, "rd_bc_validation.rds"))

checkpoint <- file.path(out, "local", "anchored_rd_checkpoint.rds")
sim <- simulate_anchored_rd(d, 500, seeds = 19:20, checkpoint_file = checkpoint, code_dir = code_dir, rd_dir = rd_dir)
resumed <- simulate_anchored_rd(d, 500, seeds = 19:20, checkpoint_file = checkpoint, code_dir = code_dir, rd_dir = rd_dir)
check(identical(sim$estimates, resumed$estimates), "completed checkpoint resumes without recalculation")
check(fails(simulate_anchored_rd(d, 500, seeds = 19:20, checkpoint_file = checkpoint,
  model_weighting = "maic", code_dir = code_dir, rd_dir = rd_dir)), "checkpoint rejects changed analysis settings")
metrics <- rd_write_simulation(sim, file.path(out, "local"))$metrics
check(nrow(metrics) == 9 && all(metrics$n_total == 2) && all(metrics$effect_scale == "RD"),
  "shared summary plus RD adapter preserve natural RD truth and status counts")
source(file.path(rd_dir, "summary/plot_rd_results.R"), encoding = "UTF-8")
plot_rd_results(file.path(out, "local", "anchored_rd.rds"), file.path(out, "plots"))
check(file.exists(file.path(out, "plots", "rd_method_comparison.png")), "RD summary exports PNG/PDF scientific comparison figures")
cluster <- rd_run_scenario(code_dir, rd_dir, file.path(out, "cluster_wald"), n = 100,
  ess_ratio = .5, seeds = 71:72, cores = 2L, ci_mode = "wald")
check(nrow(cluster$estimates) == 18 && isTRUE(cluster$complete), "two-worker cluster-style seed execution and collection complete")
check(identical(before, rd_source_hashes(code_dir, rd_dir)), "all original and RD source hashes unchanged during checks")
writeLines(c("PASS: RD checks complete", paste("Output:", out)), file.path(out, "check_results.txt"))
cat("RD CHECK OUTPUT:", out, "\n")
