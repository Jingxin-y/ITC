#!/usr/bin/env Rscript
# Standalone RD plotting: bias/SE boxplots, accuracy/coverage lollipop figures.
# Style adapted from summary/plot_cluster_results.R; base R only.
# Reads saved RD RDS files only; no source(), fitting, bootstrap or simulation.
rdplot_method_order <- c('MN_eff','RD_MAIC_sandwich','IB','RIB','RLPM','brm','brm-FC')

rdplot_metrics <- function (data, n_expected) 
{
    data$evaluation_truth <- rdplot_target_truth(data)
    allowed_point <- c("available", "boundary_truncated", "interval_failed", "se_failed", "error", "timeout", 
        "not_requested")
    keys <- unique(data[c("method", "ci_type", "component")])
    metrics <- lapply(seq_len(nrow(keys)), function(i) {
        key <- keys[i, ]
        d <- data[data$method == key$method & data$ci_type == key$ci_type & data$component == key$component & 
            data$record_complete, , drop = FALSE]
        point <- d$status %in% allowed_point & is.finite(d$estimate)
        se <- point & d$ci_type == "Wald" & is.finite(d$se) & d$se > 0 & (is.na(d$se_status) | d$se_status == 
            "available")
        ci <- point & d$status %in% c("available", "boundary_truncated") & is.finite(d$low) & is.finite(d$up) & 
            d$low <= d$up
        true <- d$evaluation_truth
        truth <- is.finite(true)
        evaluated <- point & truth
        covered <- ci & truth & d$low <= true & true <= d$up
        ncover <- sum(ci & truth)
        ncovered <- sum(covered, na.rm = TRUE)
        mean0 <- function(v) if (length(v)) 
            mean(v)
        else NA_real_
        sd0 <- function(v) if (length(v) > 1) 
            sd(v)
        else NA_real_
        sdall <- sd0(d$estimate[point])
        sdse <- sd0(d$estimate[se])
        mse <- mean0(d$se[se])
        coverage <- if (ncover) 
            ncovered/ncover
        else NA_real_
        errors <- d$estimate[evaluated] - true[evaluated]
        cbind(key, data.frame(n_expected = n_expected, n_completed_rows = nrow(d), n_point = sum(point), 
            n_se = sum(se), n_ci = sum(ci), n_boundary = sum(ci & d$status == "boundary_truncated"), 
            n_evaluated = sum(evaluated), n_coverage = ncover, n_covered = ncovered, n_not_applicable = sum(grepl("^not_applicable", 
                d$status)), mean_estimate = mean0(d$estimate[point]), mean_truth = mean0(true[evaluated]), 
            bias = mean0(errors), bias_mcse = if (length(errors) > 1) 
                sd(errors)/sqrt(length(errors))
            else NA_real_, rmse = if (length(errors)) 
                sqrt(mean(errors^2))
            else NA_real_, empirical_sd = sdall, mean_se = mse, median_se = if (any(se)) 
                median(d$se[se])
            else NA_real_, max_se = if (any(se)) 
                max(d$se[se])
            else NA_real_, empirical_sd_for_accuracy = sdse, accuracy = if (is.finite(sdse) && sdse > 
                0) 
                mse/sdse
            else NA_real_, coverage = coverage, coverage_mcse = if (ncover) 
                sqrt(coverage * (1 - coverage)/ncover)
            else NA_real_, returned_ci_fraction = sum(ci)/n_expected, returned_and_covered_fraction = if (any(truth)) 
                ncovered/n_expected
            else NA_real_, mean_width = mean0(d$up[ci] - d$low[ci]), stringsAsFactors = FALSE))
    })
    ans <- do.call(rbind, metrics)
    rownames(ans) <- NULL
    order_methods <- rdplot_method_order
    rank <- match(ans$method, order_methods)
    rank[is.na(rank)] <- length(order_methods) + 1
    ans[order(rank, ans$method, match(ans$ci_type, c("Wald", "BC"))), ]
}

rdplot_samples <- function (data, m) 
{
    data$evaluation_truth <- rdplot_target_truth(data)
    allowed_point <- c("available", "boundary_truncated", "interval_failed", "se_failed", "error", "timeout", 
        "not_requested")
    samples <- lapply(seq_len(nrow(m)), function(i) {
        d <- data[data$method == m$method[i] & data$ci_type == m$ci_type[i] & data$component == m$component[i] & 
            data$record_complete, , drop = FALSE]
        point <- d$status %in% allowed_point & is.finite(d$estimate)
        evaluated <- point & is.finite(d$evaluation_truth)
        se <- point & d$ci_type == "Wald" & is.finite(d$se) & d$se > 0 & (is.na(d$se_status) | d$se_status == 
            "available")
        list(bias = d$estimate[evaluated] - d$evaluation_truth[evaluated], se = d$se[se])
    })
    ans <- list(bias = lapply(samples, `[[`, "bias"), se = lapply(samples, `[[`, "se"))
    stopifnot(all(lengths(ans$bias) == m$n_evaluated), all(lengths(ans$se) == m$n_se))
    ans
}

rdplot_excluded_method <- function (method) 
{
    method <- toupper(gsub("_", "-", trimws(as.character(method)), fixed = TRUE))
    grepl("^BRMAD", method) | method %in% c("LP-BFC", "RLP-BFC")
}

rdplot_axis <- function (values, metric, se_max = 1, bias_max = 1) 
{
    bounds <- switch(metric, bias = c(-bias_max, bias_max), se = c(0, se_max), accuracy = c(0, 2), coverage = c(0, 
        1))
    reference <- switch(metric, bias = 0, se = NA_real_, accuracy = 1, coverage = 0.95)
    finite <- values[is.finite(values)]
    clipped <- c(any(finite < bounds[1]), any(finite > bounds[2]))
    boxes <- metric %in% c("bias", "se")
    capped <- boxes && any(clipped)
    shown <- if (capped || !boxes) 
        pmax(pmin(finite, bounds[2]), bounds[1])
    else finite
    domain <- c(shown, reference[is.finite(reference)])
    if (!length(finite) || capped) 
        limits <- bounds
    else {
        limits <- range(domain)
        if (diff(limits) == 0) 
            limits <- limits + c(-1, 1) * max(0.05, abs(limits[1]) * 0.05)
        else limits <- limits + c(-0.06, 0.06) * diff(limits)
        if (!boxes) 
            limits <- pmax(pmin(limits, bounds[2]), bounds[1])
        if (metric == "se") 
            limits[1] <- max(0, limits[1])
    }
    ticks <- pretty(limits, n = 5)
    if (!capped && length(finite)) {
        limits <- range(c(limits, ticks))
        if (!boxes) 
            limits <- pmax(pmin(limits, bounds[2]), bounds[1])
        if (metric == "se") 
            limits[1] <- max(0, limits[1])
    }
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    if (capped) 
        ticks <- sort(unique(c(ticks, bounds)))
    if (length(ticks) < 2) 
        ticks <- limits
    list(limits = limits, xlim = limits + c(-0.04, 0.04) * diff(limits), ticks = ticks, reference = reference, 
        clipped = clipped, capped = capped, bounds = bounds)
}

rdplot_metric_alpha <- function (values, metric) 
{
    target <- switch(metric, coverage = 0.95, accuracy = 1, stop("Unsupported intensity metric."))
    scale <- if (metric == "coverage") 
        0.05
    else 0.5
    distance <- abs(values - target)
    alpha <- 0.18 + 0.82 * distance/(distance + scale)
    alpha[!is.finite(distance)] <- 0.18
    alpha
}

rdplot_metric_base <- function (metric, ci_type = "Wald") 
{
    if (metric != "coverage") 
        return("#444444")
    ifelse(ci_type == "BC", "#D26920", ifelse(ci_type == "MN-ESS", "#8064A2", "#2474A6"))
}

rdplot_metric_colours <- function (values, metric, ci_type = rep("Wald", length(values))) 
{
    if (!length(values)) 
        return(character())
    base <- rep_len(rdplot_metric_base(metric, ci_type), length(values))
    alpha <- rdplot_metric_alpha(values, metric)
    vapply(seq_along(values), function(i) adjustcolor(base[i], alpha.f = alpha[i]), character(1))
}

rdplot_palette <- function (method) 
{
    cols <- c(MN_eff = "#A8A8A8", RD_MAIC_sandwich = "#7F8C8D", IB = "#D7AB87", RIB = "#83B99B", RLPM = "#ADC18B", 
        brm = "#B5A0D3", `brm-FC` = "#C19ECC")
    ans <- unname(cols[method])
    ans[is.na(ans)] <- "#AAB5C0"
    ans
}

rdplot_visible <- function (m, metric) 
{
    m <- m[!rdplot_excluded_method(m$method), , drop = FALSE]
    if (metric %in% c("bias", "se")) {
        m <- m[order(match(m$ci_type, c("Wald", "BC"))), , drop = FALSE]
        m <- m[!duplicated(m$method), , drop = FALSE]
    }
    else if (metric == "accuracy") 
        m <- m[m$ci_type == "Wald", , drop = FALSE]
    rank <- match(m$method, rdplot_method_order)
    rank[is.na(rank)] <- 99
    m[order(rank, m$method, match(m$ci_type, c("Wald", "BC"))), , drop = FALSE]
}

rdplot_draw <- function (m, metric, title, subtitle, samples, se_max = 1, bias_max = 1) 
{
    n <- nrow(m)
    x <- seq_len(n)
    boxes <- metric %in% c("bias", "se")
    conditional <- any(m$component == "conditional")
    colours <- rdplot_palette(m$method)
    values <- if (boxes) 
        unlist(samples[[metric]], use.names = FALSE)
    else m[[metric]]
    ax <- rdplot_axis(values, metric, se_max, bias_max)
    reference <- ax$reference
    yrange <- ax$limits
    span <- diff(yrange)
    ylim <- yrange + c(-0.12, 0.15) * span
    ylabels <- c(bias = "Estimation error in RD", se = "Estimated SE of RD", accuracy = if (conditional) "Mean SE / Monte Carlo SD of error" else "Mean SE / Monte Carlo SD", 
        coverage = "Coverage rate")
    titles <- c(bias = "All-method bias distributions", se = "All-method SE distributions", accuracy = "All-method SE calibration", 
        coverage = "All-method coverage rate")
    labels <- m$method
    labels[labels == "RD_MAIC_sandwich"] <- "RD-MAIC-SW"
    if (metric == "coverage") 
        labels <- paste(labels, m$ci_type, sep = " / ")
    par(mar = c(8.1, 5.2, 5.4, 1.3), mgp = c(2.7, 0.75, 0), las = 1, xpd = FALSE, family = "sans", cex = 1, 
        col.axis = "#555555", col.lab = "#202020")
    plot(NA, xlim = c(0.4, n + 0.6), ylim = ylim, xlab = "", ylab = ylabels[[metric]], xaxt = "n", yaxt = "n", 
        bty = "n", xaxs = "i", yaxs = "i")
    minor <- sort(unique(c(ax$ticks, head(ax$ticks, -1) + diff(ax$ticks)/2)))
    abline(h = minor, col = "#F0F0F0", lwd = 0.8)
    abline(h = ax$ticks, col = "#E4E4E4", lwd = 1.05)
    abline(v = x, col = "#ECECEC", lwd = 1)
    axis(2, at = ax$ticks, labels = if (metric == "coverage") 
        sprintf("%.0f%%", 100 * ax$ticks)
    else format(ax$ticks, trim = TRUE), tick = FALSE, cex.axis = 0.95)
    if (is.finite(reference)) 
        abline(h = reference, lty = 2, lwd = 1.35, col = if (metric == "bias") 
            "#3D7EAD"
        else "#333333")
    mtext(titles[[metric]], 3, line = 3.45, adj = 0, font = 2, cex = 1.65)
    mtext(title, 3, line = 2.05, adj = 0, cex = 1.03)
    mtext(subtitle, 3, line = 0.7, adj = 0, cex = 0.8, col = "#555555")
    par(xpd = NA)
    text(x, par("usr")[3] - 0.026 * diff(par("usr")[3:4]), labels, srt = 43, adj = 1, cex = 0.85, col = "#555555")
    mtext("Method", 1, line = 4.65, cex = 1.03)
    par(xpd = FALSE)
    counts <- if (metric == "bias") 
        m$n_evaluated
    else if (metric == "se") 
        m$n_se
    else if (metric == "accuracy") 
        m$n_se
    else m$n_coverage
    count_line <- paste0("n=", counts)
    fmt <- function(v) ifelse(abs(v) < 10000, sprintf("%.3f", v), formatC(v, digits = 3, format = "g"))
    if (boxes) {
        keep <- which(lengths(samples[[metric]]) > 0)
        if (length(keep)) {
            b <- boxplot(samples[[metric]][keep], plot = FALSE, range = 1.5)
            b$stats <- pmax(pmin(b$stats, ax$limits[2]), ax$limits[1])
            visible <- b$out >= ax$limits[1] & b$out <= ax$limits[2]
            b$out <- b$out[visible]
            b$group <- b$group[visible]
            bxp(b, at = x[keep], horizontal = FALSE, add = TRUE, axes = FALSE, frame.plot = FALSE, show.names = FALSE, 
                boxwex = 0.66, boxfill = adjustcolor(colours[keep], alpha.f = 0.9), boxcol = "#4D4D4D", 
                medcol = "#333333", medlwd = 2.1, whiskcol = "#555555", staplecol = "#555555", staplewex = 0, 
                whisklty = 1, outcol = adjustcolor("#444444", alpha.f = 0.42), outpch = 16, outcex = 0.35)
        }
        means <- if (metric == "bias") 
            m$bias
        else m$mean_se
        points(x, pmax(pmin(means, ax$limits[2]), ax$limits[1]), pch = 21, bg = "white", col = "#333333", 
            cex = 0.76, lwd = 0.9)
        for (side in c(-1, 1)) {
            if (metric == "se" && side < 0) 
                next
            bound <- if (side < 0) 
                -bias_max
            else if (metric == "bias") 
                bias_max
            else se_max
            clipped <- vapply(samples[[metric]], function(v) sum(side * v > side * bound), integer(1))
            ii <- which(clipped > 0)
            if (length(ii)) {
                points(x[ii], rep(bound, length(ii)), pch = 21, bg = "#D94C42", col = "#A7322A", cex = 1.9)
                text(x[ii], rep(bound, length(ii)), labels = clipped[ii], col = "white", cex = 0.62, 
                  font = 2)
            }
        }
        absent <- lengths(samples[[metric]]) == 0
        note <- "White points: full-data means. Red circles: display-clipped counts; statistics use all valid values."
    }
    else {
        ok <- is.finite(m[[metric]])
        v <- m[[metric]]
        stem_colours <- rdplot_metric_colours(v, metric, m$ci_type)
        shown <- pmax(pmin(v, ax$limits[2]), ax$limits[1])
        clipped <- ok & v > ax$limits[2]
        if (any(ok)) {
            segments(x[ok], reference, x[ok], shown[ok], lwd = 5.3, col = stem_colours[ok])
            points(x[ok], shown[ok], pch = ifelse(m$ci_type[ok] == "BC", 17, 16), col = "#222222", 
                cex = 0.92)
        }
        labs <- if (metric == "coverage") 
            sprintf("%.1f%%", 100 * v)
        else fmt(v)
        labs[clipped] <- paste0(labs[clipped], " >")
        above <- shown >= reference
        if (any(ok)) 
            text(x[ok], shown[ok] + ifelse(above[ok], 0.036, -0.04) * span, labs[ok], cex = 0.73, col = ifelse(clipped[ok], 
                "#A7322A", "#333333"))
        absent <- !ok
        note <- if (metric == "accuracy") 
            paste("Target = 1.", if (conditional) 
                "Matched replicates; denominator SD(error)."
            else "Matched SE-valid replicates; BC and MN_eff model SE not provided.", "Darker stems = larger deviation.")
        else "Target = 95%. Wald: blue; BC: orange; MN-ESS: purple. Darker stems = larger deviation. Valid returned intervals only."
    }
    if (any(absent)) 
        text(x[absent], yrange[1] + 0.06 * span, "NA", col = "#8A8A8A", cex = 0.8)
    text(x, rep(yrange[1] - 0.077 * span, n), count_line, cex = 0.64, col = "#777777")
    mtext(note, 1, line = 6.1, cex = 0.73, adj = 0)
    if (metric == "bias") 
        mtext(paste0("Bias reference = 0; critical values +/-", bias_max, if (ax$capped) 
            "; capped range activated."
        else "; adaptive axis, no observations clipped."), 1, line = 6.85, cex = 0.7, adj = 0)
    if (metric == "se") 
        mtext(paste0("SE critical value = ", se_max, if (ax$capped) 
            "; capped range activated."
        else "; adaptive axis, no observations clipped."), 1, line = 6.85, cex = 0.7, adj = 0)
}

rdplot_save <- function (metrics, folder, title, subtitle, formats, rows, se_max = 1, bias_max = 1) 
{
    items <- c("bias", "se", "accuracy", "coverage")
    draw <- function(metric) {
        mm <- rdplot_visible(metrics, metric)
        rdplot_draw(mm, metric, title, subtitle, rdplot_samples(rows, mm), se_max, bias_max)
    }
    if ("png" %in% formats) 
        for (metric in items) {
            type <- if (capabilities("cairo")) 
                "cairo"
            else if (.Platform$OS.type == "windows") 
                "windows"
            else stop("Headless PNG needs Cairo support. Use --format=pdf.")
            png(file.path(folder, paste0(metric, ".png")), width = 4200, height = 2220, res = 300, type = type)
            tryCatch(draw(metric), finally = dev.off())
        }
    if ("pdf" %in% formats) {
        pdf(file.path(folder, "metrics.pdf"), width = 14, height = 7.4, onefile = TRUE)
        tryCatch(for (metric in items) draw(metric), finally = dev.off())
    }
}

rdplot_target_truth <- function (data) 
data$evaluation_truth_saved

rdplot_usage <- function() cat(paste(
  "Usage: Rscript plot_rd_cluster_results.R RESULT_DIR [options]",
  "RESULT_DIR may be one RD run or the directory containing RD_* runs.",
  "--out=DIR          Output root; default RESULT_DIR/plots.",
  "--batch=TEXT       Select scenario folder names containing TEXT (literal match).",
  "--format=png,pdf   png, pdf, or png,pdf (default).",
  "--bias-max=NUMBER  Absolute display limit for RD errors; default 1.",
  "--se-max=NUMBER    Upper display limit for RD SE; default 1.",
  "--dry-run          List scenarios without reading results or creating plots.",
  sep = "\n"), "\n")

rdplot_options <- function(args) {
  if (!length(args) || any(args %in% c("-h", "--help"))) { rdplot_usage(); return(NULL) }
  opt <- list(input = args[1], out = NULL, batch = "", formats = c("png", "pdf"),
    bias_max = 1, se_max = 1, dry_run = FALSE)
  for (a in args[-1]) {
    if (a == "--dry-run") opt$dry_run <- TRUE
    else if (startsWith(a, "--out=")) opt$out <- substring(a, 7)
    else if (startsWith(a, "--batch=")) opt$batch <- substring(a, 9)
    else if (startsWith(a, "--format=")) opt$formats <- strsplit(substring(a, 10), ",", fixed = TRUE)[[1]]
    else if (startsWith(a, "--bias-max=")) opt$bias_max <- suppressWarnings(as.numeric(substring(a, 12)))
    else if (startsWith(a, "--se-max=")) opt$se_max <- suppressWarnings(as.numeric(substring(a, 10)))
    else stop("Unknown argument: ", a)
  }
  if (!dir.exists(opt$input)) stop("RD result directory does not exist: ", opt$input)
  opt$input <- normalizePath(opt$input, winslash = "/", mustWork = TRUE)
  if (is.null(opt$out)) opt$out <- file.path(opt$input, "plots")
  if (!nzchar(opt$out)) stop("--out must not be empty.")
  if (!length(opt$formats) || any(!opt$formats %in% c("png", "pdf"))) stop("Use --format=png,pdf, png, or pdf.")
  for (v in c("bias_max", "se_max"))
    if (length(opt[[v]]) != 1L || !is.finite(opt[[v]]) || opt[[v]] <= 0)
      stop(v, " must be a positive finite number.")
  opt
}

plot_rd_cluster_results <- function(root, out = file.path(root, "plots"), batch = "",
    formats = c("png", "pdf"), bias_max = 1, se_max = 1, dry_run = FALSE) {
  if (!dir.exists(root)) stop("RD result directory does not exist: ", root)
  runs <- if (file.exists(file.path(root, "anchored_rd.rds"))) root else
    list.dirs(root, full.names = TRUE, recursive = FALSE)
  runs <- sort(runs[grepl("^RD_", basename(runs)) | file.exists(file.path(runs, "anchored_rd.rds"))])
  if (nzchar(batch)) runs <- runs[grepl(batch, basename(runs), fixed = TRUE)]
  if (!length(runs)) {
    entries <- head(list.files(root, all.files = FALSE), 20L)
    stop("No RD scenarios matched.\nInput directory: ", normalizePath(root, winslash = "/"),
      "\nBatch filter: ", if (nzchar(batch)) batch else "<empty; all scenarios>",
      "\nTop-level entries: ", if (length(entries)) paste(entries, collapse = ", ") else "<empty or unreadable>",
      "\nUse the directory containing the original RD_* runs and anchored_rd.rds results.")
  }
  if (dry_run) { cat(paste(runs, collapse = "\n"), "\n"); return(invisible(runs)) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  index <- all_metrics <- list()
  for (run in runs) {
    input <- file.path(run, "anchored_rd.rds"); folder <- file.path(out, basename(run))
    state <- "plotted"; message <- ""; n_expected <- n_complete <- NA_integer_
    if (!file.exists(input)) {
      state <- "skipped"; message <- "No merged anchored_rd.rds result."
    } else tryCatch({
      z <- readRDS(input)
      if (!inherits(z, "anchored_rd_simulation") || !identical(z$effect_scale, "RD"))
        stop("An RD simulation result is required; RR results are rejected.")
      if (!isTRUE(z$complete)) stop("Merged result is incomplete.")
      d <- z$estimates
      if (!nrow(d) || !all(c("seed", "method", "ci_type", "component", "status", "interval_status", "se_status",
          "estimate", "se", "lower", "upper", "population_RD_AC_in_BC") %in% names(d)))
        stop("Unsupported RD result schema.")
      if (anyNA(d[c("seed", "method", "ci_type", "component", "status", "interval_status")]))
        stop("Missing seed, method or status metadata.")
      if (anyDuplicated(d[c("seed", "method", "ci_type")])) stop("Duplicate seed/method/CI rows.")
      if (any(d$component != "marginal")) stop("Marginal RD results are required.")
      d$method <- sub("-BC$", "", d$method)
      point <- d$status %in% c("available", "boundary_truncated") & is.finite(d$estimate)
      ci <- d$interval_status %in% c("available", "boundary_truncated", "BC_boundary_truncated") &
        is.finite(d$lower) & is.finite(d$upper) & d$lower <= d$upper
      # Retain valid point estimates when their interval or variance failed.
      d$status <- ifelse(point, ifelse(ci,
        ifelse(d$interval_status %in% c("boundary_truncated", "BC_boundary_truncated"),
          "boundary_truncated", "available"), "interval_failed"), "fit_failed")
      d$record_complete <- TRUE; d$low <- d$lower; d$up <- d$upper
      d$evaluation_truth_saved <- d$population_RD_AC_in_BC
      n_complete <- length(unique(d$seed)); n_expected <- length(z$requested_seeds)
      spec_file <- file.path(run, "submission.txt")
      if (file.exists(spec_file)) {
        spec <- readLines(spec_file, warn = FALSE); value <- spec[startsWith(spec, "nsim=")]
        if (length(value) == 1) n_expected <- as.integer(sub("^nsim=", "", value))
      }
      if (!is.finite(n_expected) || n_expected < n_complete) stop("Invalid requested repeat count.")
      m <- rdplot_metrics(d, n_expected)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      title <- paste("BC marginal A-C RD |", z$design$event, "| n =", z$n_AC,
        "| ESS ratio =", format(z$design$ess_ratio, digits = 3))
      a <- z$settings$analysis
      adjustment <- if (is.null(a$estimation_mode)) "scores" else a$estimation_mode
      brm_model <- if (adjustment == "none") "intercept_only" else
        if (is.null(a$brm_effect_model)) "saturated" else a$brm_effect_model
      subtitle <- paste("Adjustment:", adjustment, "| fit weights:", a$model_weighting,
        "| BRM:", brm_model,
        "| source:", z$design$ess_source, "| EM:", z$design$em_level,
        "| calibration:", paste(a$weight_covariates, collapse = ","),
        "| complete", n_complete, "/", n_expected)
      rdplot_save(m, folder, title, subtitle, formats, d, se_max, bias_max)
      m$effect_scale <- "RD"; m$run_name <- basename(run); m$brm_effect_model <- brm_model
      write.csv(m, file.path(folder, "metrics.csv"), row.names = FALSE)
      writeLines(c("Target: BC population marginal A-C RD, on the natural RD scale.",
        "Bias boxplots use individual estimate-minus-truth errors, one per method per seed.",
        "Accuracy = mean(SE) / SD(estimates), calculated on the SAME SE-valid Wald replicates.",
        "BC and MN_eff provide no model SE and are excluded from the accuracy figure.",
        "Coverage uses valid returned intervals. Also inspect returned_ci_fraction and returned_and_covered_fraction.",
        "Display clipping affects pictures only, never the metrics. Different runs are not pooled."),
        file.path(folder, "analysis_notes.txt"))
      all_metrics[[length(all_metrics) + 1L]] <- m
    }, error = function(e) { state <<- "failed"; message <<- conditionMessage(e) })
    index[[length(index) + 1L]] <- data.frame(run_name = basename(run), output_dir = folder,
      status = state, n_expected = n_expected, n_complete = n_complete, message = message)
    cat(toupper(state), ": ", folder, if (nzchar(message)) paste0(" | ", message) else "", "\n", sep = "")
  }
  ix <- do.call(rbind, index)
  write.csv(ix, file.path(out, "rd_plot_run_index.csv"), row.names = FALSE)
  if (length(all_metrics)) write.csv(do.call(rbind, all_metrics),
    file.path(out, "rd_plot_metrics_all_runs.csv"), row.names = FALSE)
  if (any(ix$status == "failed")) stop("Some RD scenarios failed; see rd_plot_run_index.csv.")
  if (!any(ix$status == "plotted")) stop("No merged RD results were plotted.")
  invisible(list(index = ix, metrics = all_metrics))
}

if (sys.nframe() == 0L) {
  if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
  opt <- rdplot_options(commandArgs(TRUE))
  if (!is.null(opt)) plot_rd_cluster_results(opt$input, opt$out, opt$batch,
    opt$formats, opt$bias_max, opt$se_max, opt$dry_run)
}
