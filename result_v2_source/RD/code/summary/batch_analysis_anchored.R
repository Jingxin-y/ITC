## Rscript summary/batch_analysis_anchored.R input.rds-or-directory [output_dir] [truth.csv]
## 旧 CSV 绘图入口完整保存在 summary/legacy/batch_analysis_anchored.R。
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript summary/batch_analysis_anchored.R input.rds-or-directory [output_dir] [truth.csv]")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
source(file.path(dirname(script), "summarize_joint_simulation.R"), encoding = "UTF-8")
input <- if (dir.exists(args[1])) file.path(args[1], "anchored_joint.rds") else args[1]
result <- readRDS(input)
truth <- if (length(args) >= 3L) read.csv(args[3], stringsAsFactors = FALSE) else result$analysis_truth
summary <- summarize_anchored_simulation(result, truth)
out_dir <- if (length(args) >= 2L) args[2] else dirname(input)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary$metrics, file.path(out_dir, "anchored_joint_summary.csv"), row.names = FALSE)
write.csv(summary$method_status, file.path(out_dir, "anchored_joint_status.csv"), row.names = FALSE)
write.csv(summary$parameter_status, file.path(out_dir, "anchored_joint_parameter_status.csv"), row.names = FALSE)
if (!is.null(result$nonconverged_estimates))
  write.csv(result$nonconverged_estimates,file.path(out_dir,"anchored_joint_nonconverged.csv"),row.names=FALSE)
cat("Saved summaries to", normalizePath(out_dir), "\n")
