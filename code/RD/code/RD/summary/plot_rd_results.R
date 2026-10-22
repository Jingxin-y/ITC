plot_rd_results <- function(input, output_dir = file.path(dirname(input), "plots")) {
  if (dir.exists(input)) input <- file.path(input, "anchored_rd.rds")
  z <- readRDS(input)
  if (!inherits(z, "anchored_rd_simulation") || !identical(z$effect_scale, "RD"))
    stop("An RD simulation RDS is required; RR results are rejected.")
  m <- summarize_anchored_rd_simulation(z)$metrics
  if (!nrow(m)) stop("No RD results to plot.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(m, file.path(output_dir, "rd_plot_metrics.csv"), row.names = FALSE)
  metrics <- c("bias", "rmse", "coverage", "estimate_failure_rate", "empirical_sd", "mean_se")
  labels <- c("Bias (RD)", "RMSE (RD)", "95% interval coverage", "Estimate failure rate", "Empirical SD (RD)", "Mean model SE (RD)")
  draw <- function() {
    old <- par(mfrow = c(2, 3), mar = c(8, 5.5, 2.5, 1), mgp = c(4, 1, 0), oma = c(0, 0, 3, 0))
    on.exit(par(old))
    colors <- ifelse(grepl("-BC$", m$method), "#c07929", "#247c92")
    for (j in seq_along(metrics)) {
      values <- m[[metrics[j]]]; limits <- range(c(0, values[is.finite(values)]))
      if (metrics[j] %in% c("coverage", "estimate_failure_rate")) limits <- c(0, 1)
      if (diff(limits) == 0) limits <- limits + c(-.01, .01)
      centers <- barplot(ifelse(is.finite(values), values, 0), names.arg = m$method,
        col = colors, las = 2, cex.names = .72, border = NA, ylim = limits,
        ylab = labels[j], main = labels[j])
      if (any(!is.finite(values))) text(centers[!is.finite(values)], rep(limits[2] * .06, sum(!is.finite(values))), "NA", cex = .7)
      abline(h = if (metrics[j] == "coverage") .95 else 0, lty = 2, col = "#888888")
    }
    title <- paste("BC population marginal A-C RD |", z$design$event, z$design$ess_source,
      z$design$em_level, "| n =", z$n_AC, "ESS =", format(z$design$ess_ratio, digits = 3),
      "| repeats =", length(z$seeds))
    mtext(title, outer = TRUE, cex = 1.1)
  }
  png(file.path(output_dir, "rd_method_comparison.png"), width = 2100, height = 1400, res = 180)
  tryCatch(draw(), finally = dev.off())
  pdf(file.path(output_dir, "rd_method_comparison.pdf"), width = 12, height = 8)
  tryCatch(draw(), finally = dev.off())
  writeLines(c("Scale: natural RD; target: BC population marginal A-C effect.",
    "Coverage uses available verified intervals; failures are reported separately.",
    "MN_eff SE is not a model SE: CI-width diagnostic is saved in the result table.",
    "BC is conditional parametric-bootstrap inversion, not a finite-sample exact guarantee.",
    "Wald intervals and LPM counterfactual predictions are not clipped."), file.path(output_dir, "README.txt"))
  invisible(m)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Usage: Rscript RD/summary/plot_rd_results.R RD_RESULT_DIR [OUTPUT_DIR] [SHARED_CODE_DIR]")
  own <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  rd_dir <- dirname(dirname(normalizePath(own, winslash = "/")))
  if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
  source(file.path(rd_dir, "load_anchored_rd.R"), encoding = "UTF-8")
  source_anchored_rd(if (length(args) >= 3) args[3] else dirname(rd_dir), rd_dir)
  plot_rd_results(args[1], if (length(args) >= 2) args[2] else file.path(args[1], "plots"))
}
