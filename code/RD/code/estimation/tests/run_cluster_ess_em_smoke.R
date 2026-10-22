#!/usr/bin/env Rscript
## Real pipeline smoke test; reduced bootstrap budgets, not a coverage study.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript run_cluster_ess_em_smoke.R CODE_DIR OUTPUT_DIR [CORES]")
code_dir <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
run_dir <- args[2L]
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
if (length(cores) != 1L || is.na(cores) || cores < 1L) stop("Invalid core count.")
if (.Platform$OS.type != "unix" && cores > 1L) stop("Use one core on Windows.")
allocation <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (!is.na(allocation) && cores > allocation) stop("Core count exceeds the Slurm allocation.")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
source_anchored(code_dir)
source(file.path(code_dir, "estimation/cluster_all_methods.R"), encoding = "UTF-8")
source(file.path(code_dir, "estimation/cluster_grid.R"), encoding = "UTF-8")

cases <- data.frame(
  event = rep(c("common", "rare"), each = 3L),
  ess_source = rep(c("diffuse", "em_driven", "nonem_driven"), 2L),
  em_level = rep(c("low", "high", "high"), 2L),
  weight_mode = rep(c("all", "em", "nonem"), 2L),
  n = 200L, ess_ratio = c(.25, .15, .15, .5, .5, .5), seed = 1L,
  stringsAsFactors = FALSE
)
cases$case <- with(cases, paste(event, ess_source, em_level, weight_mode,
                               paste0("n", n), paste0("ess", 100 * ess_ratio), sep = "_"))
write.csv(cases, file.path(run_dir, "smoke_cases.csv"), row.names = FALSE)

run_case <- function(i) {
  cell <- cases[i, ]
  Sys.setenv(MAIC_DGP = "ess_em", MAIC_EVENT = cell$event, MAIC_ESS_SOURCE = cell$ess_source,
             MAIC_EM_LEVEL = cell$em_level, MAIC_WEIGHT_MODE = cell$weight_mode,
             MAIC_MODEL_WEIGHTING = "unit", MAIC_EST_MODE = "scores",
             MAIC_EST_EM = "z1,z2,z3", MAIC_EST_NON_EM = "z4,z5,z6")
  folder <- file.path(run_dir, "scenarios", cell$case)
  dir.create(file.path(folder, "logs"), recursive = TRUE, showWarnings = FALSE)
  fields <- c(nsim = 1L, first_seed = cell$seed, ci_mode = "smoke", event = cell$event,
              estimation_mode = "scores", n = cell$n, ess_ratio = cell$ess_ratio, cores = 1L,
              dgp_mode = "ess_em", ess_source = cell$ess_source, em_level = cell$em_level,
              weight_mode = cell$weight_mode, model_weighting = "unit")
  writeLines(paste(names(fields), fields, sep = "="), file.path(folder, "submission.txt"))
  cat("SMOKE start:", cell$case, "\n"); flush.console()
  maic_run_scenario(code_dir, folder, cell$n, cell$ess_ratio, "smoke", cores = 1L)
  maic_collect(folder)
  checkpoint <- readRDS(file.path(folder, "by_seed", sprintf("seed_%06d.rds", cell$seed)))
  d <- checkpoint$result$simulation$design
  expected_cols <- paste0("z", switch(cell$weight_mode, all = 1:6, em = 1:3, nonem = 4:6))
  stopifnot(checkpoint$complete, is.null(checkpoint$error), checkpoint$settings$version == 10L,
    checkpoint$settings$counterfactual_risk_check == "finite_positive_standardized_means",
            identical(checkpoint$settings$dgm, "webster2024_rr_ess_em_v1"),
            identical(checkpoint$settings$ess_source, cell$ess_source),
            identical(checkpoint$settings$em_level, cell$em_level),
            identical(checkpoint$settings$weight_covariates, expected_cols),
            identical(checkpoint$settings$model_weighting, "unit"),
            all(checkpoint$result$model_weights$weight == 1),
            identical(checkpoint$result$variables$standardization_weighting, "maic"),
            isTRUE(d$calibration$matches_target_BC),
            abs(d$ess_ratio - cell$ess_ratio) < 1e-10)
  exact <- lapply(c("brm", "brm-FC"), function(method) {
    ci <- checkpoint$exact[[method]]
    if (is.null(ci) || is.null(ci$diagnostics$counts) ||
        ci$diagnostics$counts$profile.fits < 1L || ci$diagnostics$counts$bootstrap.draws < 1L)
      stop(cell$case, ": ", method, " did not execute real bootstrap/profile calculations.")
    data.frame(case = cell$case, event = cell$event, method = method,
      status = ci$status, interval_available = isTRUE(ci$interval.available),
      low = ci$low, up = ci$up,
      profile_fits = ci$diagnostics$counts$profile.fits,
      bootstrap_draws = ci$diagnostics$counts$bootstrap.draws,
      stringsAsFactors = FALSE)
  })
  combined <- do.call(rbind, exact)
  cat("SMOKE checked:", cell$case, "\n"); flush.console()
  combined
}

outcomes <- parallel::mclapply(seq_len(nrow(cases)), function(i) {
  tryCatch(list(rows = run_case(i), error = NULL),
           error = function(e) list(rows = NULL, error = conditionMessage(e)))
}, mc.cores = min(cores, nrow(cases)), mc.preschedule = FALSE, mc.set.seed = FALSE)
bad <- vapply(outcomes, function(x) !is.list(x) || !is.null(x$error), logical(1L))
if (any(bad)) {
  failures <- data.frame(case = cases$case[bad], message = vapply(outcomes[bad],
    function(x) if (is.list(x)) x$error else "Worker failed without a result.", character(1L)))
  write.csv(failures, file.path(run_dir, "smoke_failures.csv"), row.names = FALSE)
  print(failures)
  stop("Smoke test failed; inspect smoke_failures.csv and the per-case checkpoints.")
}
rows <- do.call(rbind, lapply(outcomes, `[[`, "rows"))
write.csv(rows, file.path(run_dir, "smoke_exact_checks.csv"), row.names = FALSE)
print(rows)
for (event in c("common", "rare")) {
  if (!any(rows$event == event & rows$interval_available &
           rows$status %in% c("available", "boundary_truncated")))
    stop("No returned BC interval in ", event, "; inspect smoke_exact_checks.csv.")
}
message <- c("PASS: six complete data-generation/estimation/checkpoint/collection cases.",
             "PASS: common/rare, three ESS sources, low/high EM and all/em/nonem calibration.",
             "PASS: unit fitting weights and MAIC marginal standardization.",
             "PASS: real BRM and BRM-FC bootstrap/profile CI calculations; returned intervals in both events.",
             "Reduced bootstrap counts are for execution testing only; they do not establish coverage.")
writeLines(message, file.path(run_dir, "SMOKE_PASS.txt"))
cat(paste(message, collapse = "\n"), "\n")
