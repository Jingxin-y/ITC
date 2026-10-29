args <- commandArgs(trailingOnly = TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
code_dir <- if (length(args)) args[1] else dirname(rd_dir)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8"); source_anchored_rd(code_dir, rd_dir)
fails <- function(x) inherits(tryCatch(x, error = identity), "error")
control <- list(brm_max_step = 5L, brm_optim_maxit = 10L, fc_max_step = 10L, glm_maxit = 20L)
sim <- generate_webster_rd_scenario(100, event = "rare", seed = 29)
sim$data.IPD$data$y <- 0
z <- get.compare.anchored.rd(sim$data.IPD, sim$data.AD, fit_control = control)
stopifnot(nrow(z$estimates) == 9,
  all(z$estimates$status[z$estimates$method %in% c("RLPM", "MN_eff", "RD_MAIC_sandwich")] == "available"),
  all(z$estimates$estimate[z$estimates$method %in% c("RLPM", "MN_eff", "RD_MAIC_sandwich")] == 0))
cat("PASS: zero-event outcomes retain available LPM/MN/weighted-RD results and separate failed fits\n")
sim <- generate_webster_rd_scenario(200, seed = 33)
sim$data.IPD$data$z1 <- 1; sim$data.IPD$data$z2 <- 0; sim$data.IPD$data$z3 <- 1
z <- get.compare.anchored.rd(sim$data.IPD, sim$data.AD, weight_covariates = character(0),
  score_standardization = "fixed", fit_control = control)
stopifnot(all(z$estimates$status[z$estimates$method %in% c("brm", "brm-FC")] == "fit_failed"),
  all(z$estimates$status[z$estimates$method %in% c("MN_eff", "RD_MAIC_sandwich")] == "available"))
cat("PASS: absent EM levels fail only methods requiring an identifiable effect design\n")
sim$data.AD$anchored.mean.x["z1"] <- 2
z <- get.compare.anchored.rd(sim$data.IPD, sim$data.AD, score_standardization = "fixed")
stopifnot(all(z$estimates$status == "weighting_failed"), nrow(z$estimates) == 9)
cat("PASS: invalid MAIC target preserves all nine failure rows\n")
sim <- generate_webster_rd_scenario(500, seed = 19)
stopifnot(fails(get.compare.anchored.rd(sim$data.IPD, sim$data.AD, param = "RR")),
  fails(get.compare.anchored.rd(sim$data.IPD, sim$data.AD, bc_control = list(unknown = 1))),
  fails(get.compare.anchored.rd(sim$data.IPD, sim$data.AD, bc_control = list(pars = 1))))
cat("PASS: wrong scales and unsupported/reserved BC controls are rejected\n")
z <- get.compare.anchored.rd(sim$data.IPD, sim$data.AD, run_exact = TRUE,
  seed = 19, bc_control = list(time.limit.sec = 1e-12))
idx <- z$estimates$ci_type == "BC"
stopifnot(all(z$estimates$status[idx] == "available"), all(z$estimates$interval_status[idx] == "timeout"),
  all(is.finite(z$estimates$estimate[idx])))
cat("PASS: BC timeout retains independent point estimates\n")

## Test the constrained likelihood gradient independently of the RD map derivatives.
dat <- sim$data.IPD$data; va <- cbind(1, as.numeric(dat$S_E == 1), as.numeric(dat$S_E == 2), as.numeric(dat$S_E == 3))
vb <- va; pars <- c(.1, .05, .1, .2, -.5, .1, .2, .3); w <- exp(.4 * dat$z1)
problem <- .rd_exact_problem(dat$y, dat$x, va, vb, w, pars, va, vb, w,
  intercept = 1L, alpha.bound = 15, beta.bound = 15, optim.maxit = 100,
  optim.reltol = 1e-10, score.tol = 1e-4, root.tol = 1e-10, root.maxit = 200,
  constraint.tol = 1e-7, eps = 1e-12, check.deadline = function(stage) NULL)
psi <- .12; analytic <- problem$evaluate(problem$start, psi, dat$y)
ng <- vapply(seq_along(problem$start), function(j) {
  up <- down <- problem$start; up[j] <- up[j] + 1e-5; down[j] <- down[j] - 1e-5
  (problem$evaluate(up, psi, dat$y)$value - problem$evaluate(down, psi, dat$y)$value) / 2e-5
}, numeric(1))
stopifnot(max(abs(ng - analytic$gradient)) < 1e-4, abs(analytic$constraint.error) < 1e-7)
cat("PASS: constrained saturated-RD profile gradient agrees with numerical differentiation\n")
rare <- run.anchored.rd(webster_rd_scenario_design("rare", .15, "em_driven", "high"),
  100, seed = 7, fit_control = control)
stopifnot(nrow(rare$estimates) == 9,
  all(rare$estimates$status != "pending"), all(rare$estimates$effect_scale == "RD"))
cat("PASS: rare/low-ESS small sample exports complete independent method statuses\n")
