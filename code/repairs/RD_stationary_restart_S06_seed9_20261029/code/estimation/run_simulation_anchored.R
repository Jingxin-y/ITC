## Rscript estimation/run_simulation_anchored.R config.R [output_dir]
## config.R 定义 simulation_args，以及可选 truth；不再使用旧 event/hypothesis 参数。
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript estimation/run_simulation_anchored.R config.R [output_dir]")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
code_dir <- dirname(dirname(script))
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
source_anchored(code_dir)
config <- new.env(parent = globalenv())
source(args[1], local = config, encoding = "UTF-8")
if (!exists("simulation_args", config, inherits = FALSE))
  stop("Configuration must define simulation_args for simulate_anchored().")
out_dir <- if (length(args) >= 2L) args[2] else "results_joint"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
simulation_args <- config$simulation_args
if (is.null(simulation_args$checkpoint_file))
  simulation_args$checkpoint_file <- file.path(out_dir, "anchored_joint.rds")
result <- do.call(simulate_anchored, simulation_args)
truth <- if (exists("truth", config, inherits = FALSE)) config$truth else NULL
if (is.null(truth)) truth <- result$reference_truth
result$analysis_truth <- truth
saveRDS(result, file.path(out_dir, "anchored_joint.rds"))
write.csv(result$estimates, file.path(out_dir, "anchored_joint_estimates.csv"), row.names = FALSE)
write.csv(result$estimates[result$estimates$component=="marginal",,drop=FALSE],
          file.path(out_dir,"anchored_joint_marginal_exact.csv"),row.names=FALSE)
write.csv(result$nonconverged_estimates,file.path(out_dir,"anchored_joint_nonconverged.csv"),row.names=FALSE)
saveRDS(list(estimates=result$nonconverged_estimates,fits=result$nonconverged_fits,
  design=result$design,n_AC=result$n_AC,n_BC=result$n_BC,fit_control=result$fit_control),
  file.path(out_dir,"anchored_joint_nonconverged.rds"))
summary <- summarize_anchored_simulation(result, truth)
write.csv(summary$metrics, file.path(out_dir, "anchored_joint_summary.csv"), row.names = FALSE)
write.csv(summary$method_status, file.path(out_dir, "anchored_joint_status.csv"), row.names = FALSE)
write.csv(summary$parameter_status, file.path(out_dir, "anchored_joint_parameter_status.csv"), row.names = FALSE)
cat("Saved", length(result$seeds), "replicates to", normalizePath(out_dir), "\n")
