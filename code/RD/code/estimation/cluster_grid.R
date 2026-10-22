## Grid orchestration only: estimation and CI algorithms live in the existing helpers.
maic_grid_scenarios <- function(run_dir) {
  tab <- read.csv(file.path(run_dir, "scenarios.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(tab) == 9L, !anyDuplicated(tab$scenario),
    all(grepl("^n(100|200|500)_ess(015|050|080)$", tab$scenario)),
    setequal(tab$n, c(100L, 200L, 500L)), setequal(tab$ess_ratio, c(.15, .5, .8)),
    nrow(unique(tab[c("n", "ess_ratio")])) == 9L,
    all(tab$scenario == sprintf("n%d_ess%03d", tab$n, round(100 * tab$ess_ratio))))
  tab
}

maic_run_grid_task <- function(code_dir, run_dir, task_id, ci_mode = c("full", "smoke")) {
  ci_mode <- match.arg(ci_mode)
  if (!all(c("n", "ess_ratio") %in% names(formals(maic_run_seed))))
    stop("Upload the updated estimation/cluster_all_methods.R with n and ess_ratio arguments.")
  scenarios <- maic_grid_scenarios(run_dir)
  tasks <- read.csv(file.path(run_dir, "tasks.csv"), stringsAsFactors = FALSE)
  task <- tasks[tasks$task_id == task_id, , drop = FALSE]
  stopifnot(nrow(task) == 1L, task$scenario %in% scenarios$scenario,
    task$first_seed >= 1L, task$last_seed >= task$first_seed)
  scenario <- scenarios[match(task$scenario, scenarios$scenario), , drop = FALSE]
  stopifnot(task$n == scenario$n, task$ess_ratio == scenario$ess_ratio)
  folder <- file.path(run_dir, "scenarios", task$scenario)
  cat("Task", task_id, "scenario", task$scenario, "seeds", task$first_seed, "to", task$last_seed, "\n")
  errors <- data.frame(seed = integer(), message = character())
  for (seed in seq.int(task$first_seed, task$last_seed)) {
    error <- tryCatch({
      maic_run_seed(code_dir, folder, seed, ci_mode, n = task$n, ess_ratio = task$ess_ratio)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(error)) {
      errors <- rbind(errors, data.frame(seed = seed, message = error))
      cat("Seed", seed, "failed:", error, "\n")
      write.csv(errors, file.path(run_dir, "logs", sprintf("task_%04d_failures.csv", task_id)), row.names = FALSE)
    }
  }
  if (nrow(errors)) stop(nrow(errors), " replicate(s) failed in task ", task_id, "; other replicates were processed.")
  invisible(NULL)
}

maic_collect_grid <- function(run_dir) {
  scenarios <- maic_grid_scenarios(run_dir)
  rows <- coefficients <- missing <- incomplete <- progress <- list()
  add_scenario <- function(d, name) { d$scenario <- rep(name, nrow(d)); d }
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios$scenario[i]
    folder <- file.path(run_dir, "scenarios", scenario)
    cat("Collecting", scenario, "\n")
    completed <- maic_collect(folder)
    absent <- read.csv(file.path(folder, "missing_seeds.csv"))
    partial <- read.csv(file.path(folder, "incomplete_seeds.csv"))
    missing[[i]] <- add_scenario(absent, scenario)
    incomplete[[i]] <- add_scenario(partial, scenario)
    file <- file.path(folder, "combined_results.csv")
    recorded <- 0L
    if (file.exists(file)) {
      d <- read.csv(file, stringsAsFactors = FALSE)
      if (any(d$n_AC != scenarios$n[i] | d$n_BC != scenarios$n[i] |
              abs(d$ess_ratio - scenarios$ess_ratio[i]) > 1e-12))
        stop("Result settings do not match scenario ", scenario)
      rows[[i]] <- add_scenario(d, scenario)
      recorded <- length(unique(d$seed))
    }
    cf <- file.path(folder, "combined_coefficients.csv")
    if (file.exists(cf)) coefficients[[i]] <- add_scenario(read.csv(cf, stringsAsFactors = FALSE), scenario)
    progress[[i]] <- data.frame(scenario = scenario, n = scenarios$n[i], ess_ratio = scenarios$ess_ratio[i],
      recorded = recorded, complete = completed, incomplete = nrow(partial), missing = nrow(absent))
  }
  bind <- function(x) do.call(rbind, Filter(Negate(is.null), x))
  combined <- bind(rows)
  if (!is.null(combined)) {
    write.csv(combined, file.path(run_dir, "combined_results.csv"), row.names = FALSE)
    counts <- aggregate(list(n_rows = rep(1L, nrow(combined))),
      combined[c("scenario", "n_AC", "n_BC", "ess_ratio", "method", "ci_type", "status")], sum)
    write.csv(counts, file.path(run_dir, "status_counts.csv"), row.names = FALSE)
  }
  cf <- bind(coefficients)
  if (!is.null(cf)) write.csv(cf, file.path(run_dir, "combined_coefficients.csv"), row.names = FALSE)
  write.csv(bind(missing), file.path(run_dir, "missing_seeds.csv"), row.names = FALSE)
  write.csv(bind(incomplete), file.path(run_dir, "incomplete_seeds.csv"), row.names = FALSE)
  summary <- bind(progress)
  write.csv(summary, file.path(run_dir, "grid_progress.csv"), row.names = FALSE)
  print(summary)
  invisible(summary)
}

# A single Slurm job owns exactly one scenario and forks one process per seed.
# Each process writes only its own seed files, so checkpoints do not collide.
maic_run_scenario <- function(code_dir, run_dir, n, ess_ratio,
                              ci_mode = c("full", "smoke"), cores = 104L) {
  ci_mode <- match.arg(ci_mode)
  if (.Platform$OS.type != "unix" && cores > 1L)
    stop("Parallel scenario jobs require Linux/Unix.")
  if (!all(c("n", "ess_ratio") %in% names(formals(maic_run_seed))))
    stop("Upload the updated estimation/cluster_all_methods.R with n and ess_ratio arguments.")
  spec <- readLines(file.path(run_dir, "submission.txt"), warn = FALSE)
  field <- function(key) {
    values <- sub(paste0("^", key, "="), "", spec[startsWith(spec, paste0(key, "="))])
    if (length(values) != 1L) stop("Missing or duplicate submission field: ", key)
    values
  }
  if (!identical(as.integer(field("n")), n) ||
      !isTRUE(all.equal(as.numeric(field("ess_ratio")), ess_ratio)) ||
      !identical(field("ci_mode"), ci_mode) ||
      !identical(as.integer(field("cores")), cores))
    stop("Worker arguments do not match scenario submission settings.")
  if (!identical(field("event"), tolower(Sys.getenv("MAIC_EVENT", "common"))))
    stop("Worker event does not match scenario submission settings.")
  if (!identical(field("estimation_mode"), tolower(Sys.getenv("MAIC_EST_MODE", "scores"))))
    stop("Worker estimation mode does not match scenario submission settings.")
  optional_env <- c(dgp_mode = "MAIC_DGP", ess_source = "MAIC_ESS_SOURCE",
                    em_level = "MAIC_EM_LEVEL", weight_mode = "MAIC_WEIGHT_MODE",
                    model_weighting = "MAIC_MODEL_WEIGHTING")
  defaults <- c(dgp_mode = "legacy", ess_source = "diffuse", em_level = "high",
                weight_mode = "all", model_weighting = "unit")
  for (key in names(optional_env)) if (any(startsWith(spec, paste0(key, "=")))) {
    if (!identical(field(key), tolower(Sys.getenv(optional_env[key], defaults[key]))))
      stop("Worker ", key, " does not match scenario submission settings.")
  }
  if (length(cores) != 1L || is.na(cores) || cores < 1L || cores != floor(cores))
    stop("cores must be a positive integer.")
  allocation <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (!is.na(allocation) && cores > allocation)
    stop("Requested R workers exceed SLURM_CPUS_PER_TASK.")
  nsim <- as.integer(field("nsim"))
  first_seed <- as.integer(field("first_seed"))
  if (is.na(nsim) || nsim < 1L || is.na(first_seed) || first_seed < 1L ||
      first_seed > .Machine$integer.max - 1000000L - nsim)
    stop("Invalid seed range in submission.txt.")
  seeds <- seq.int(first_seed, length.out = nsim)
  cat("Scenario n", n, "ESS", ess_ratio, "seeds", first_seed, "to", tail(seeds, 1L),
      "on", cores, "cores\n")
  outcomes <- parallel::mclapply(seeds, function(seed) {
    error <- tryCatch({
      maic_run_seed(code_dir, run_dir, seed, ci_mode, n = n, ess_ratio = ess_ratio)
      NULL
    }, error = function(e) conditionMessage(e))
    list(seed = seed, error = error)
  }, mc.cores = min(cores, nsim), mc.preschedule = FALSE, mc.set.seed = FALSE)
  bad <- vapply(outcomes, function(x) !is.list(x) || !is.null(x$error), logical(1L))
  if (any(bad)) {
    errors <- data.frame(seed = seeds[bad], message = vapply(outcomes[bad],
      function(x) if (is.list(x)) x$error else "Worker failed without a result", character(1L)))
    write.csv(errors, file.path(run_dir, "logs", "seed_failures.csv"), row.names = FALSE)
    stop(nrow(errors), " replicate(s) failed; other seeds were processed.")
  }
  cat("All", nsim, "replicates processed.\n")
  invisible(NULL)
}
