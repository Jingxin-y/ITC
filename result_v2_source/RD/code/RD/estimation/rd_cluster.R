rd_smoke_bc_control <- function() list(nsim.coarse = 2L, nsim.confirm = 2L,
  nsim.refine = 2L, nsim.p.value = 2L, thres.dicho = .05, bisection.max.step = 8L,
  optim.maxit = 100L, optim.retry.maxit = 300L, time.limit.sec = 60)
rd_run_seed <- function(code_dir, rd_dir, run_dir, seed, n, ess_ratio,
    event = "common", ess_source = "diffuse", em_level = "high", weight_mode = "all",
    model_weighting = "unit", ci_mode = "smoke", bc_control = NULL,
    estimation_mode = "scores", brm_effect_model = "linear_atanh") {
  if (!weight_mode %in% c("all", "em", "nonem") || !ci_mode %in% c("smoke", "full", "wald"))
    stop("Invalid weight/CI mode.")
  covs <- switch(weight_mode, all = paste0("z", 1:6), em = paste0("z", 1:3), nonem = paste0("z", 4:6))
  if (is.null(bc_control)) bc_control <- if (ci_mode == "smoke") rd_smoke_bc_control() else list()
  design <- webster_rd_scenario_design(event, ess_ratio, ess_source, em_level)
  path <- file.path(run_dir, "by_seed", sprintf("rd_seed_%08d.rds", seed))
  result <- simulate_anchored_rd(design, n, n, seeds = seed, run_exact = ci_mode != "wald",
    weight_covariates = covs, model_weighting = model_weighting, estimation_mode = estimation_mode,
    brm_effect_model = brm_effect_model,
    checkpoint_file = path, bc_control = bc_control, code_dir = code_dir, rd_dir = rd_dir)
  cat("RD seed", seed, "complete;", sum(result$estimates$interval_status == "available"), "intervals available\n")
  invisible(path)
}
rd_run_scenario <- function(code_dir, rd_dir, run_dir, n = 100L, ess_ratio = .5,
    seeds = 1:2, cores = 1L, event = "common", ess_source = "diffuse", em_level = "high",
    weight_mode = "all", model_weighting = "unit", ci_mode = "smoke",
    estimation_mode = "scores", brm_effect_model = "linear_atanh") {
  estimation_mode <- match.arg(estimation_mode, c("scores", "none"))
  brm_effect_model <- match.arg(brm_effect_model, c("linear_atanh", "saturated"))
  rd_validate_sizes(n, n); rd_validate_seeds(seeds)
  if (!is.numeric(cores) || length(cores) != 1 || !is.finite(cores) || cores < 1 || cores != floor(cores))
    stop("cores must be a positive integer.")
  args <- list(code_dir = code_dir, rd_dir = rd_dir, run_dir = run_dir, n = n, ess_ratio = ess_ratio,
    event = event, ess_source = ess_source, em_level = em_level, weight_mode = weight_mode,
    model_weighting = model_weighting, ci_mode = ci_mode, estimation_mode = estimation_mode,
    brm_effect_model = brm_effect_model)
  job <- function(seed, args) do.call(rd_run_seed, c(args, list(seed = seed)))
  if (cores == 1 || length(seeds) == 1) lapply(seeds, job, args = args) else {
    cl <- parallel::makeCluster(min(cores, length(seeds)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(code_dir, rd_dir) {
      if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
      source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
      source_anchored_rd(code_dir, rd_dir); NULL
    }, code_dir, rd_dir)
    parallel::parLapply(cl, seeds, job, args = args)
  }
  rd_collect(run_dir)
}
rd_collect <- function(run_dir) {
  paths <- sort(list.files(file.path(run_dir, "by_seed"), "^rd_seed_.*\\.rds$", full.names = TRUE))
  if (!length(paths)) stop("No RD seed checkpoints found.")
  records <- lapply(paths, readRDS); first <- records[[1]]$settings; first$requested_seeds <- NULL
  for (z in records) {
    current <- z$settings; current$requested_seeds <- NULL
    if (!inherits(z, "anchored_rd_simulation") || !identical(current, first))
      stop("Mixed effect scale, source hashes, DGM or analysis settings in RD collection.")
    if (!isTRUE(z$complete) || length(z$seeds) != 1) stop("Incomplete or non-single-seed checkpoint in collection.")
  }
  seeds <- unlist(lapply(records, function(z) z$seeds))
  if (anyDuplicated(seeds)) stop("Duplicate seed checkpoints.")
  combined <- records[[1]]
  combined$estimates <- do.call(rbind, lapply(records, function(z) z$estimates))
  combined$coefficient_estimates <- do.call(rbind, lapply(records, function(z) z$coefficient_estimates))
  combined$reference_truth <- do.call(rbind, lapply(records, function(z) z$reference_truth))
  combined$diagnostics <- unlist(lapply(records, function(z) z$diagnostics), recursive = FALSE)
  combined$seeds <- combined$requested_seeds <- seeds
  combined$settings$requested_seeds <- seeds
  combined$parameter_schema <- unique(combined$estimates[c("method", "component", "parameter", "term")])
  rd_write_simulation(combined, run_dir)
  invisible(combined)
}
