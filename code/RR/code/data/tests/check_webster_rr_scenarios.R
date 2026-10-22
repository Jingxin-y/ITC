## From anchored_v1: Rscript data/tests/check_webster_rr_scenarios.R
code_dir <- Sys.getenv("MAIC_CODE_DIR", unset = ".")
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
source_anchored_data(code_dir)
close_to <- function(a, b, tol = 1e-10)
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
bad <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
z <- as.matrix(expand.grid(rep(list(0:1), 6)))
E <- (1 - z[, 1]) + z[, 2] + (1 - z[, 3])
P <- z[, 4] + (1 - z[, 5]) + z[, 6]
index <- 1L + E + 4L * P
count <- 0L
expected_targets <- rbind(
  common = c(p_C = .418596226792393, p_A = .514669087904676, p_B = .60),
  rare = c(p_C = .0488309131792698, p_A = .10, p_B = .15))

for (event in c("common", "rare")) for (level in c("low", "high")) {
  base <- webster_rr_scenario_design(event = event, em_level = level)
  close_to(webster_rr_scenario_targets(event), expected_targets[event, ])
  close_to(base$parameters$baseline_rr_E, .59)
  close_to(base$parameters$baseline_rr_P, 1.08)
  stopifnot(base$scenario$version == 2L,
            identical(base$scenario$risk_target_version, "BC_CAB_20261015"))
  lambda <- if (level == "low") .15 else .60
  close_to(base$coefficient_truth[, "lambda"], c(lambda, lambda / 2))
  close_to(base$scenario$lambda_AB, lambda / 2)
  for (src in c("diffuse", "em_driven", "nonem_driven")) for (r in c(.15, .25, .5, .8, 1)) {
    d <- webster_rr_scenario_design(event, r, src, level)
    stopifnot(identical(d$dgm, "webster2024_rr_ess_em_v2"),
              identical(d$parameters, base$parameters),
              identical(d$risk_table, base$risk_table),
              isTRUE(d$calibration$matches_target_BC),
              identical(d$calibration$target_population, "BC"))
    close_to(d$ess_ratio, r)
    expected <- switch(src, diffuse = c(r, sqrt(r), sqrt(r)),
                       em_driven = c(r, r, 1), nonem_driven = c(r, 1, r))
    close_to(d$relative_ess, expected)
    if (src == "em_driven") close_to(d$q[4:6], d$e[4:6])
    if (src == "nonem_driven") close_to(d$q[1:3], d$e[1:3])
    if (r == 1) close_to(d$q, d$e)

    ## Independent integration over all 64 raw covariate configurations.
    tab <- d$risk_table[index, ]
    stopifnot(all(as.matrix(tab[c("p_A", "p_B", "p_C")]) > 0),
              all(as.matrix(tab[c("p_A", "p_B", "p_C")]) < 1))
    for (population in c("AC", "BC")) {
      prob <- if (population == "AC") d$q else d$e
      mass <- apply(sweep(z, 2, prob, function(x, p) p^x * (1 - p)^(1 - x)), 1, prod)
      risks <- colSums(tab[c("p_C", "p_A", "p_B")] * mass)
      close_to(risks, d$marginal_truth[population, names(risks)])
      if (population == "BC") {
        close_to(risks, expected_targets[event, names(risks)])
        close_to(d$marginal_truth[population, "log_RR_AB"],
                 log(expected_targets[event, "p_A"] / expected_targets[event, "p_B"]))
      }
    }
    if (src == "nonem_driven") close_to(d$marginal_truth["AC", c("log_RR_AC", "log_RR_BC", "log_RR_AB")],
                                       d$marginal_truth["BC", c("log_RR_AC", "log_RR_BC", "log_RR_AB")])
    count <- count + 1L
  }
}
cat("PASS:", count, "designs; exact ESS decomposition, legal 64-cell risks, fixed BC means and marginal truths\n")

a <- generate_webster_rr_scenario(200, 240, "rare", .15, "em_driven", "high", seed = 27)
b <- generate_webster_rr_scenario(200, 240, "rare", .15, "em_driven", "high", seed = 27)
stopifnot(identical(a, b), identical(names(a), c("data.IPD", "data.AD.or", "data.AD", "design")))
for (name in c("data.IPD", "data.AD.or")) {
  trial <- a[[name]]
  n <- if (name == "data.IPD") 200 else 240
  close_to(trial$count[c("Na0", "Na1")], c(n / 2, n / 2))
  stopifnot(nrow(trial$data) == n, all(trial$data$y %in% 0:1))
  tab <- a$design$risk_table[trial$truth$risk_table_row, ]
  active <- if (name == "data.IPD") tab$p_A else tab$p_B
  close_to(trial$truth$observed_risk, ifelse(trial$data$x == 1, active, tab$p_C))
}
close_to(a$data.AD$anchored.mean.x, colMeans(a$data.AD.or$data[paste0("z", 1:6)]))
cat("PASS: reproducible data, balanced trial arms, shared control risks and existing AD interface\n")

for (r in list(0, -1, 1.1, NA_real_, Inf, "0.5")) bad(webster_rr_scenario_design(ess_ratio = r))
bad(webster_rr_scenario_design(ess_source = "wrong"))
bad(webster_rr_scenario_design(em_level = "wrong"))
bad(webster_rr_scenario_design(event = "wrong"))
for (n in list(0, 3, NA_real_, Inf, c(100, 200))) bad(generate_webster_rr_scenario(n))
for (seed in list(-1, 1.5, NA_real_, Inf)) bad(generate_webster_rr_scenario(100, seed = seed))
cat("PASS: invalid settings rejected\n")
