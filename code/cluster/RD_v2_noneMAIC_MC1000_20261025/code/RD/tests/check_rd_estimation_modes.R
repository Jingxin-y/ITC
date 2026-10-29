args <- commandArgs(trailingOnly = TRUE)
own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
code_dir <- if (length(args)) args[1] else dirname(rd_dir)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
source_anchored_rd(code_dir, rd_dir)
check <- function(x, message) { if (!isTRUE(x)) stop(message); cat("PASS:", message, "\n") }
fails <- function(x) inherits(tryCatch(x, error = identity), "error")
design <- webster_rd_scenario_design("common", .5, "diffuse", "high")
scores <- run.anchored.rd(design, 500, seed = 19)
none <- run.anchored.rd(design, 500, seed = 19, estimation_mode = "none")
raw <- c(paste0("z", 1:6), "x", "y")
check(identical(scores$data$data.IPD$data[raw], none$data$data.IPD$data[raw]) &&
      identical(scores$data$data.AD$anchored.mean.x[paste0("z", 1:6)],
                none$data$data.AD$anchored.mean.x[paste0("z", 1:6)]) &&
      identical(scores$data$data.AD$anchored.cov.x[paste0("z", 1:6), paste0("z", 1:6)],
                none$data$data.AD$anchored.cov.x[paste0("z", 1:6), paste0("z", 1:6)]) &&
      identical(scores$weights$weight, none$weights$weight),
      "scores/none preserve generated raw data, BC summaries and all-six-covariate MAIC weights")
check(identical(none$variables$weight_covariates, paste0("z", 1:6)) &&
      identical(none$variables$estimation_mode, "none") &&
      !length(none$variables$model_covariates) && !length(none$variables$effect_modifiers) &&
      !length(none$variables$brm_effect_modifiers) &&
      all(vapply(none$fits[c("brm", "brm-FC")], function(f)
        isTRUE(f$converged) && ncol(f$va) == 1L && ncol(f$vb) == 1L, logical(1))),
      "none uses intercept-only BRM/FC effect and OP, while keeping full MAIC calibration")
dat <- none$data$data.IPD$data
observed_rd <- mean(dat$y[dat$x == 1]) - mean(dat$y[dat$x == 0])
oracle_fit <- run.anchored.rd(design, 500, seed = 19, estimation_mode = "none",
  fit_control = list(brm_optim_reltol = 1e-14, brm_thres = 1e-14))
check(ncol(none$fits$IB$design$X) == 2L && ncol(none$fits$RLPM$design$X) == 2L &&
      all(abs(none$estimates$estimate[none$estimates$method %in% c("IB", "RIB", "RLPM")] -
        observed_rd) < 1e-10) &&
      abs(oracle_fit$estimates$estimate[oracle_fit$estimates$method == "brm"] - observed_rd) < 1e-6,
      "unit-weight none regressions and tightly optimized BRM match independent arm risk difference")
methods <- c("MN_eff", "RD_MAIC_sandwich")
cols <- c("method", "estimate", "se", "lower", "upper", "p.value")
check(identical(scores$estimates[scores$estimates$method %in% methods, cols],
                none$estimates[none$estimates$method %in% methods, cols]),
      "MAIC arm-mean methods are unchanged by estimation mode")
check(fails(run.anchored.rd(design, 500, seed = 19, estimation_mode = "unknown")),
      "unknown estimation modes are rejected")
out <- Sys.getenv("RD_TEST_OUTPUT", tempfile("rd_modes_"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(out, "none_checkpoint.rds")
sim <- simulate_anchored_rd(design, 500, seeds = 19L, estimation_mode = "none",
  checkpoint_file = checkpoint, code_dir = code_dir, rd_dir = rd_dir)
resumed <- simulate_anchored_rd(design, 500, seeds = 19L, estimation_mode = "none",
  checkpoint_file = checkpoint, code_dir = code_dir, rd_dir = rd_dir)
check(identical(sim$estimates, resumed$estimates) &&
      identical(sim$settings$analysis$estimation_mode, "none") &&
      all(sim$estimates$estimation_mode == "none") &&
      fails(simulate_anchored_rd(design, 500, seeds = 19L, estimation_mode = "scores",
        checkpoint_file = checkpoint, code_dir = code_dir, rd_dir = rd_dir)),
      "checkpoint records estimation mode, resumes none and rejects scores/none mixing")
parallel_result <- rd_run_scenario(code_dir, rd_dir, file.path(out, "parallel_none"),
  n = 500, ess_ratio = .5, seeds = 19:20, cores = 2L, ci_mode = "wald", estimation_mode = "none")
check(nrow(parallel_result$estimates) == 18L && isTRUE(parallel_result$complete) &&
      all(parallel_result$estimates$estimation_mode == "none") &&
      identical(parallel_result$settings$analysis$estimation_mode, "none") &&
      all(vapply(parallel_result$diagnostics, function(z)
        identical(z$brm_effect_basis, "intercept only") && identical(z$brm_score_levels, 0), logical(1))),
      "two-worker execution forwards none mode and retains intercept-only diagnostics")
smoke <- rd_run_scenario(code_dir, rd_dir, file.path(out, "smoke_none"),
  n = 500, ess_ratio = .5, seeds = 19L, cores = 1L, ci_mode = "smoke", estimation_mode = "none")
check(all(smoke$estimates$interval_status == "available") &&
      all(smoke$estimates$p_status == "available") &&
      identical(smoke$settings$analysis$estimation_mode, "none"),
      "actual none smoke includes real BC inversion and nine available intervals/P values")
writeLines("PASS: scores/none estimation-mode checks complete", file.path(out, "check_rd_estimation_modes.txt"))
