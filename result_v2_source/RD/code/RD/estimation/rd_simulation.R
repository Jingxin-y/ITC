rd_source_hashes <- function(code_dir = dirname(.anchored_rd_module_dir), rd_dir = .anchored_rd_module_dir) {
  parent <- c("load_anchored.R", unlist(lapply(c("data", "estimation", "summary"), function(d)
    file.path(d, list.files(file.path(code_dir, d), pattern = "\\.[Rr]$", recursive = TRUE)))))
  own <- list.files(rd_dir, pattern = "\\.(R|sh)$", recursive = TRUE)
  paths <- c(file.path(code_dir, parent), file.path(rd_dir, own))
  setNames(unname(tools::md5sum(paths)), c(paste0("shared/", parent), paste0("RD/", own)))
}
rd_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid()); saveRDS(object, tmp)
  if (.Platform$OS.type == "windows" && file.exists(path) && !file.remove(path))
    stop("Cannot replace checkpoint: ", path)
  if (!file.rename(tmp, path)) stop("Cannot finalize checkpoint: ", path)
}
rd_simulation_settings <- function(design, n_AC, n_BC, seeds, analysis,
    code_dir = dirname(.anchored_rd_module_dir), rd_dir = .anchored_rd_module_dir) {
  list(version = 1L, effect_scale = "RD", target = "BC marginal A-C RD", methods = rd_method_plan(),
    probability_mapping = .anchored_rd_solver_env$probability_mapping,
    dgm = design$dgm, design = design, n_AC = n_AC, n_BC = n_BC,
    requested_seeds = seeds, analysis = analysis, source_md5 = rd_source_hashes(code_dir, rd_dir))
}
simulate_anchored_rd <- function(design, n_AC, n_BC = n_AC, seeds = seq_len(100),
    run_exact = FALSE, checkpoint_file = NULL, resume = TRUE,
    model_weighting = c("unit", "maic"), weight_covariates = paste0("z", 1:6),
    score_spec = NULL, score_standardization = c("bc_estimated", "fixed"),
    model_covariates = NULL, effect_modifiers = NULL, fit_control = list(),
    bc_control = list(), exact_seed_offset = 1000000L,
    code_dir = dirname(.anchored_rd_module_dir), rd_dir = .anchored_rd_module_dir,
    estimation_mode = c("scores", "none"),
    brm_effect_model = c("linear_atanh", "saturated")) {
  rd_validate_sizes(n_AC, n_BC); rd_validate_seeds(seeds); rd_validate_seeds(seeds + exact_seed_offset)
  if (!identical(design$effect_scale, "RD")) stop("RD design required.")
  estimation_mode <- match.arg(estimation_mode)
  brm_effect_model <- match.arg(brm_effect_model)
  model_weighting <- match.arg(model_weighting); score_standardization <- match.arg(score_standardization)
  validate_rd_bc_control(bc_control)
  analysis <- list(run_exact = run_exact, model_weighting = model_weighting,
    weight_covariates = weight_covariates, score_spec = score_spec,
    score_standardization = score_standardization, model_covariates = model_covariates,
    effect_modifiers = effect_modifiers, fit_control = anchored_fit_control(fit_control),
    bc_control = bc_control, exact_seed_offset = exact_seed_offset,
    estimation_mode = estimation_mode, brm_effect_model = brm_effect_model)
  settings <- rd_simulation_settings(design, n_AC, n_BC, seeds, analysis, code_dir, rd_dir)
  rows <- coefficients <- diagnostics <- vector("list", length(seeds)); completed <- 0L
  snapshot <- function(k) {
    dat <- if (k) do.call(rbind, rows[seq_len(k)]) else data.frame()
    truth <- if (k) unique(dat[c("r", "method", "component", "term", "population_RD_AC_in_BC")]) else NULL
    if (!is.null(truth)) names(truth)[names(truth) == "population_RD_AC_in_BC"] <- "truth"
    structure(list(settings = settings, estimates = dat,
      coefficient_estimates = if (k) do.call(rbind, coefficients[seq_len(k)]) else data.frame(),
      diagnostics = if (k) diagnostics[seq_len(k)] else list(),
      parameter_schema = if (k) unique(dat[c("method", "component", "parameter", "term")]) else NULL,
      reference_truth = truth, design = design, n_AC = n_AC, n_BC = n_BC,
      seeds = if (k) seeds[seq_len(k)] else numeric(0), requested_seeds = seeds,
      complete = k == length(seeds), effect_scale = "RD", variables = analysis),
      class = c("anchored_rd_simulation", "anchored_simulation", "list"))
  }
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    if (!resume) stop("Checkpoint already exists; choose another output or resume=TRUE.")
    old <- readRDS(checkpoint_file)
    if (!inherits(old, "anchored_rd_simulation") || !identical(old$settings, settings))
      stop("Checkpoint settings, effect scale or source hashes differ.")
    completed <- length(old$seeds)
    if (!identical(old$seeds, seeds[seq_len(completed)])) stop("Checkpoint seeds are not a valid prefix.")
    if (completed) for (i in seq_len(completed)) {
      rows[[i]] <- old$estimates[old$estimates$r == seeds[i], , drop = FALSE]
      coefficients[[i]] <- old$coefficient_estimates[old$coefficient_estimates$r == seeds[i], , drop = FALSE]
      diagnostics[[i]] <- old$diagnostics[[i]]
    }
  }
  if (completed == length(seeds)) return(snapshot(completed))
  if (!is.null(checkpoint_file) && !file.exists(checkpoint_file)) rd_atomic_save(snapshot(0L), checkpoint_file)
  for (i in seq.int(completed + 1L, length(seeds))) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(do.call(run.anchored.rd, c(list(design = design, n_AC = n_AC, n_BC = n_BC,
      seed = seeds[i]), analysis)), error = function(e) {
      if (!is.null(checkpoint_file)) {
        partial <- snapshot(i - 1L); partial$interrupted_at_seed <- seeds[i]; partial$error <- conditionMessage(e)
        rd_atomic_save(partial, checkpoint_file)
      }
      stop(e)
    })
    d <- fit$estimates; d$r <- d$seed <- seeds[i]; d$event <- design$event
    d$em_level <- design$em_level; d$ess_source <- design$ess_source
    d$ess_ratio <- design$ess_ratio; d$n_AC <- n_AC; d$n_BC <- n_BC; d$dgm <- design$dgm
    d$model_weighting <- model_weighting; d$estimation_mode <- estimation_mode
    d$sample_ESS <- if (all(is.finite(fit$weights$weight))) sum(fit$weights$weight)^2 / sum(fit$weights$weight^2) else NA_real_
    rows[[i]] <- d
    cf <- fit$coefficient_estimates
    if (nrow(cf)) { cf$r <- cf$seed <- seeds[i] }
    coefficients[[i]] <- cf
    diagnostics[[i]] <- list(seed = seeds[i], warnings = fit$warnings, errors = fit$errors,
      weight_diagnostics = fit$weights$diagnostics, score_spec = fit$variables$score_spec,
      estimation_mode = fit$variables$estimation_mode,
      brm_effect_model = fit$variables$brm_effect_model,
      brm_effect_basis = fit$variables$brm_effect_basis, brm_score_levels = fit$variables$brm_score_levels,
      bc = fit$bc, fit_convergence = lapply(fit$fits[c("brm", "brm-FC")], function(z)
        if (is.null(z)) NULL else list(converged = z$converged, optimization = z$optimization,
          fc = z$fc.diagnostics, last_optimizer = z$last_optimizer)),
      elapsed_seconds = proc.time()[["elapsed"]] - t0)
    if (!is.null(checkpoint_file)) rd_atomic_save(snapshot(i), checkpoint_file)
  }
  snapshot(length(seeds))
}
summarize_anchored_rd_simulation <- function(result) {
  if (!inherits(result, "anchored_rd_simulation") || !identical(result$effect_scale, "RD"))
    stop("RD simulation required.")
  out <- summarize_anchored_simulation(result)
  m <- out$metrics; d <- result$estimates
  if (nrow(m)) {
    m$rmse <- vapply(seq_len(nrow(m)), function(i) {
      z <- d[d$method == m$method[i] & d$component == m$component[i] & d$term == m$term[i], , drop = FALSE]
      ok <- z$status == "available" & is.finite(z$estimate) & is.finite(z$population_RD_AC_in_BC)
      if (any(ok)) sqrt(mean((z$estimate[ok] - z$population_RD_AC_in_BC[ok])^2)) else NA_real_
    }, numeric(1))
    m$interval_failure_rate <- 1 - m$n_interval / m$n_total
    m$coverage_all_replicates <- ifelse(m$n_coverage > 0, m$coverage * m$n_coverage / m$n_total, NA_real_)
    m$effect_scale <- "RD"; out$metrics <- m
  }
  out
}
rd_write_simulation <- function(result, output_dir) {
  if (!inherits(result, "anchored_rd_simulation")) stop("RD simulation required.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rd_atomic_save(result, file.path(output_dir, "anchored_rd.rds"))
  write.csv(result$estimates, file.path(output_dir, "anchored_rd_estimates.csv"), row.names = FALSE)
  write.csv(result$coefficient_estimates, file.path(output_dir, "anchored_rd_coefficients.csv"), row.names = FALSE)
  summary <- summarize_anchored_rd_simulation(result)
  write.csv(summary$metrics, file.path(output_dir, "anchored_rd_summary.csv"), row.names = FALSE)
  write.csv(summary$method_status, file.path(output_dir, "anchored_rd_status.csv"), row.names = FALSE)
  write.csv(summary$parameter_status, file.path(output_dir, "anchored_rd_parameter_status.csv"), row.names = FALSE)
  invisible(summary)
}
