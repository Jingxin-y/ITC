rm(list = ls())
gc()

raw_args <- commandArgs(trailingOnly = TRUE)
positional_args <- raw_args[!grepl("^--", raw_args)]
flag_args <- raw_args[grepl("^--", raw_args)]

input_dir <- if (length(positional_args) >= 1) positional_args[[1]] else "D:/UofT/code/RRRDOR/ITC/result/anchored"
plot_dir <- if (length(positional_args) >= 2) positional_args[[2]] else file.path(input_dir, "plots")
filters <- list(
  param = if (length(positional_args) >= 3 && positional_args[[3]] != "all") positional_args[[3]] else NA_character_,
  event = if (length(positional_args) >= 4 && positional_args[[4]] != "all") positional_args[[4]] else NA_character_,
  hypothesis = if (length(positional_args) >= 5 && positional_args[[5]] != "all") positional_args[[5]] else NA_character_,
  n = if (length(positional_args) >= 6 && positional_args[[6]] != "all") as.integer(positional_args[[6]]) else NA_integer_
)
allow_main_rlp_fallback <- "--allow-main-rlp-fallback" %in% flag_args

local_lib <- Sys.getenv("ANCHOR_R_LIB", "C:/tmp/Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required for plotting. ",
    "Install it first, for example: install.packages('ggplot2')"
  )
}

required_R <- c(200, 400, 600, 800, 1000)
method_order <- c(
  "brm", "CMH.an", "RCMH", "lb.an", "lp.an", "rlp.an",
  "GC", "brm.an", "brm.ad.an", "brm.bc.an", "brm.adbc.an"
)
data_names <- c(
  "data.brm", "data.CMH.an", "data.RCMH", "data.lb.an", "data.lp.an",
  "data.rlp.an", "data.GC", "data.brm.an", "data.brm.ad.an",
  "data.exact.an", "data.exact.ad.an"
)
names(data_names) <- method_order

est_result <- function(df, para_true) {
  est <- df$estimate
  se <- df$se
  low <- df$low
  up <- df$up

  bias <- mean(est, na.rm = TRUE) - para_true
  se_est <- mean(se, na.rm = TRUE) / sqrt(sum(!is.na(se)))
  sd_est <- mean(se, na.rm = TRUE)
  sd_mc <- sd(est, na.rm = TRUE)
  acc <- sd_est / sd_mc
  cov <- mean((low < para_true) & (up > para_true), na.rm = TRUE)
  p <- mean(df$p <= 0.05, na.rm = TRUE)

  c(bias = bias, se = se_est, acc = acc, coverage = cov, p = p)
}

alpha_true_for <- function(param, event, hypothesis) {
  if (param != "RR") {
    stop(sprintf("Unsupported param: %s", param))
  }

  if (event == "common") {
    if (hypothesis == "alternative") {
      return(0.4)
    }
    return(0)
  }

  if (event == "rare") {
    if (hypothesis == "alternative") {
      return(0.7)
    }
    return(0)
  }

  stop(sprintf("Unsupported event: %s", event))
}

parse_result_files <- function(input_dir) {
  files <- list.files(
    input_dir,
    pattern = "^ITC_simulation_results_anchored_.*\\.csv$",
    full.names = FALSE
  )
  files <- files[!grepl("_LBP\\.csv$", files)]

  pattern <- "^ITC_simulation_results_anchored_([^_]+)_([^_]+)_([^_]+)_n_([0-9]+)_R_([0-9]+)_ess_([0-9.]+)\\.csv$"
  parsed <- regexec(pattern, files)
  pieces <- regmatches(files, parsed)
  pieces <- pieces[lengths(pieces) == 7]

  if (length(pieces) == 0) {
    return(data.frame())
  }

  do.call(rbind, lapply(pieces, function(x) {
    data.frame(
      file = x[[1]],
      param = x[[2]],
      event = x[[3]],
      hypothesis = x[[4]],
      n = as.integer(x[[5]]),
      R = as.integer(x[[6]]),
      ess = as.numeric(x[[7]]),
      stringsAsFactors = FALSE
    )
  }))
}

find_complete_scenarios <- function(file_index) {
  if (nrow(file_index) == 0) {
    return(data.frame())
  }

  keys <- unique(file_index[c("param", "event", "hypothesis", "n", "ess")])
  complete <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    matches <- file_index[
      file_index$param == key$param &
        file_index$event == key$event &
        file_index$hypothesis == key$hypothesis &
        file_index$n == key$n &
        file_index$ess == key$ess,
    ]
    missing_R <- setdiff(required_R, unique(matches$R))
    rlp_file <- rlp_file_path(key)
    missing_RLP <- !file.exists(rlp_file)
    data.frame(
      key,
      complete = length(missing_R) == 0 && (!missing_RLP || allow_main_rlp_fallback),
      missing_R = paste(missing_R, collapse = ","),
      missing_RLP = missing_RLP,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, complete)
  out[order(out$param, out$event, out$hypothesis, out$n, out$ess), ]
}

filter_scenarios <- function(scenarios) {
  if (nrow(scenarios) == 0) {
    return(scenarios)
  }

  keep <- rep(TRUE, nrow(scenarios))
  if (!is.na(filters$param)) {
    keep <- keep & scenarios$param == filters$param
  }
  if (!is.na(filters$event)) {
    keep <- keep & scenarios$event == filters$event
  }
  if (!is.na(filters$hypothesis)) {
    keep <- keep & scenarios$hypothesis == filters$hypothesis
  }
  if (!is.na(filters$n)) {
    keep <- keep & scenarios$n == filters$n
  }
  scenarios[keep, , drop = FALSE]
}

rlp_file_path <- function(scenario) {
  file.path(
    input_dir,
    sprintf(
      "ITC_simulation_results_anchored_RLP_%s_%s_%s_n_%d_R_1000_ess_%s.csv",
      scenario$param, scenario$event, scenario$hypothesis, scenario$n,
      format_ess(scenario$ess)
    )
  )
}

split_methods <- function(data, scenario) {
  out <- list(
    data.brm = data[(seq_len(nrow(data)) %% 11) == 1, ],
    data.brm.an = data[(seq_len(nrow(data)) %% 11) == 2, ],
    data.lb.an = data[(seq_len(nrow(data)) %% 11) == 3, ],
    data.lp.an = data[(seq_len(nrow(data)) %% 11) == 4, ],
    data.rlp.an = data[(seq_len(nrow(data)) %% 11) == 5, ],
    data.CMH.an = data[(seq_len(nrow(data)) %% 11) == 6, ],
    data.RCMH = data[(seq_len(nrow(data)) %% 11) == 7, ],
    data.brm.ad.an = data[(seq_len(nrow(data)) %% 11) == 8, ],
    data.GC = data[(seq_len(nrow(data)) %% 11) == 9, ],
    data.exact.an = data[(seq_len(nrow(data)) %% 11) == 10, ],
    data.exact.ad.an = data[(seq_len(nrow(data)) %% 11) == 0, ]
  )

  used_lbp <- FALSE
  if (scenario$n == 50) {
    lbp_file <- file.path(
      input_dir,
      sprintf(
        "ITC_simulation_results_anchored_%s_%s_%s_n_%d_R_1000_ess_%s_LBP.csv",
        scenario$param, scenario$event, scenario$hypothesis, scenario$n,
        format_ess(scenario$ess)
      )
    )
    if (file.exists(lbp_file)) {
      data_lbp <- read.csv(lbp_file)
      out$data.lb.an <- data_lbp[(seq_len(nrow(data_lbp)) %% 11) == 3, ]
      out$data.lp.an <- data_lbp[(seq_len(nrow(data_lbp)) %% 11) == 4, ]
      used_lbp <- TRUE
    }
  }

  rlp_file <- rlp_file_path(scenario)
  used_rlp_file <- FALSE
  if (file.exists(rlp_file)) {
    data_rlp <- read.csv(rlp_file)
    out$data.rlp.an <- data_rlp[(seq_len(nrow(data_rlp)) %% 11) == 5, ]
    used_rlp_file <- TRUE
  } else if (!allow_main_rlp_fallback) {
    stop(sprintf("Missing RLP file: %s", rlp_file))
  }

  out$used_lbp <- used_lbp
  out$used_rlp_file <- used_rlp_file
  out
}

finite_rows <- function(df) {
  numeric_cols <- vapply(df, is.numeric, logical(1))
  df[apply(df[, numeric_cols, drop = FALSE], 1, function(x) all(is.finite(x))), ]
}

format_ess <- function(x) {
  sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

scenario_id <- function(scenario) {
  sprintf(
    "%s_%s_%s_n_%d_ess_%s",
    scenario$param, scenario$event, scenario$hypothesis, scenario$n,
    format_ess(scenario$ess)
  )
}

matrix_to_long <- function(x, value_name) {
  out <- data.frame(
    method = factor(rep(colnames(x), each = nrow(x)), levels = method_order),
    value = as.vector(x),
    stringsAsFactors = FALSE
  )
  names(out)[2] <- value_name
  out
}

clip_limits <- function(x, lower = -Inf, upper = Inf, include = NULL) {
  finite <- x[is.finite(x)]
  if (length(finite) == 0) {
    return(NULL)
  }

  clip_low <- is.finite(lower) && any(finite < lower)
  clip_high <- is.finite(upper) && any(finite > upper)
  if (!clip_low && !clip_high) {
    return(NULL)
  }

  visible <- finite
  if (clip_low) {
    visible <- visible[visible >= lower]
  }
  if (clip_high) {
    visible <- visible[visible <= upper]
  }
  visible <- c(visible, include)
  visible <- visible[is.finite(visible)]

  low <- if (clip_low) lower else if (length(visible) > 0) min(visible) else NA_real_
  high <- if (clip_high) upper else if (length(visible) > 0) max(visible) else NA_real_

  if (!is.finite(low) || !is.finite(high) || low >= high) {
    if (is.finite(lower) && is.finite(upper) && lower < upper) {
      return(c(lower, upper))
    }
    center <- finite[1]
    return(center + c(-0.5, 0.5))
  }

  pad <- diff(c(low, high)) * 0.06
  if (!clip_low) {
    low <- low - pad
  }
  if (!clip_high) {
    high <- high + pad
  }
  c(low, high)
}

method_counts <- function(x, predicate) {
  counts <- sapply(seq_len(ncol(x)), function(i) sum(predicate(x[, i]), na.rm = TRUE))
  data.frame(
    method = factor(colnames(x), levels = method_order),
    count = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

plot_palette <- function(n) {
  grDevices::hcl.colors(n, "Dark 3")
}

plot_estimates <- function(est, alpha_true, scenario, output_file) {
  lower <- -4
  upper <- 4
  est_long <- matrix_to_long(est, "estimate")
  y_limits <- clip_limits(est_long$estimate, lower = lower, upper = upper, include = alpha_true)

  p <- ggplot2::ggplot(est_long, ggplot2::aes(x = method, y = estimate, fill = method)) +
    ggplot2::geom_boxplot(color = "grey20", width = 0.8, outlier.size = 0.5, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = plot_palette(length(method_order))) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::geom_hline(yintercept = alpha_true, linetype = "dashed", color = "steelblue") +
    ggplot2::stat_summary(
      fun = mean, geom = "point", shape = 21, size = 1.5,
      fill = "white", color = "black", na.rm = TRUE
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 16),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = sprintf("Monte Carlo distributions of estimated log(RR)s\n%s", scenario_id(scenario)),
      x = "Method",
      y = "Estimated log(RR)"
    )

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    y_gap <- diff(y_limits) * 0.04

    hi <- method_counts(est, function(z) !is.na(z) & z > upper)
    hi <- hi[hi$count > 0, ]
    if (nrow(hi) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = hi, ggplot2::aes(x = method, y = upper, size = count),
          inherit.aes = FALSE, shape = 21, fill = "red", color = "red", alpha = 0.55
        ) +
        ggplot2::geom_text(
          data = hi, ggplot2::aes(x = method, y = upper - y_gap, label = count),
          inherit.aes = FALSE, size = 4, color = "blue"
        )
    }

    lo <- method_counts(est, function(z) !is.na(z) & z < lower)
    lo <- lo[lo$count > 0, ]
    if (nrow(lo) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = lo, ggplot2::aes(x = method, y = lower, size = count),
          inherit.aes = FALSE, shape = 21, fill = "red", color = "red", alpha = 0.55
        ) +
        ggplot2::geom_text(
          data = lo, ggplot2::aes(x = method, y = lower + y_gap, label = count),
          inherit.aes = FALSE, size = 4, color = "blue"
        )
    }

    p <- p + ggplot2::scale_size_area(max_size = 20)
  }

  ggplot2::ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}

plot_se <- function(se, scenario, output_file) {
  upper <- 4
  se_long <- matrix_to_long(se, "SE")
  se_plot <- se_long[is.finite(se_long$SE), ]
  y_limits <- clip_limits(se_plot$SE, lower = 0, upper = upper)

  p <- ggplot2::ggplot(se_plot, ggplot2::aes(x = method, y = SE, fill = method)) +
    ggplot2::geom_boxplot(color = "grey20", width = 0.8, outlier.size = 0.5, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = plot_palette(length(method_order))) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::stat_summary(
      fun = mean, geom = "point", shape = 21, size = 1.5,
      fill = "white", color = "black", na.rm = TRUE
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 16),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = sprintf("Monte Carlo distributions of standard errors\n%s", scenario_id(scenario)),
      x = "Method",
      y = "Estimated SE"
    )

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    y_gap <- diff(y_limits) * 0.04
    removed <- method_counts(se, function(z) !is.finite(z) | z > upper)
    removed <- removed[removed$count > 0, ]
    if (nrow(removed) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = removed, ggplot2::aes(x = method, y = upper, size = count),
          inherit.aes = FALSE, shape = 21, fill = "red", color = "red", alpha = 0.55
        ) +
        ggplot2::geom_text(
          data = removed, ggplot2::aes(x = method, y = upper - y_gap, label = count),
          inherit.aes = FALSE, size = 4, color = "blue"
        ) +
        ggplot2::scale_size_area(max_size = 20)
    }
  }

  ggplot2::ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}

plot_metric <- function(result, metric, target, scenario, output_file, lower = -Inf, upper = Inf) {
  y <- as.numeric(result[metric, method_order])
  metric_df <- data.frame(
    method = factor(method_order, levels = method_order),
    value = y,
    dev = abs(y - target),
    stringsAsFactors = FALSE
  )
  y_limits <- clip_limits(metric_df$value, lower = lower, upper = upper, include = target)
  base_range <- if (is.null(y_limits)) range(c(metric_df$value, target), finite = TRUE) else y_limits
  y_gap <- diff(base_range) * 0.04
  metric_df$label_y <- ifelse(metric_df$value >= target, metric_df$value + y_gap, metric_df$value - y_gap)
  if (!is.null(y_limits)) {
    metric_df$label_y <- pmin(pmax(metric_df$label_y, y_limits[1] + y_gap), y_limits[2] - y_gap)
  }

  p <- ggplot2::ggplot(metric_df, ggplot2::aes(x = method, y = value)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = method, y = target, yend = value, color = dev),
      linewidth = 3, lineend = "round"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = sprintf("%.3f", value), hjust = ifelse(value >= target, 0, 1)),
      size = 4
    ) +
    ggplot2::geom_hline(yintercept = target, linetype = "dashed") +
    ggplot2::scale_x_discrete(limits = method_order, drop = FALSE) +
    ggplot2::scale_color_gradient(low = "grey80", high = "red4") +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 16),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = sprintf("%s by method\n%s", metric, scenario_id(scenario)),
      x = "Method",
      y = metric
    )

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  ggplot2::ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}

analyze_scenario <- function(scenario) {
  files <- file.path(
    input_dir,
    sprintf(
      "ITC_simulation_results_anchored_%s_%s_%s_n_%d_R_%d_ess_%s.csv",
      scenario$param, scenario$event, scenario$hypothesis, scenario$n,
      required_R, format_ess(scenario$ess)
    )
  )
  data <- do.call(rbind, lapply(files, read.csv))
  methods <- split_methods(data, scenario)

  methods$data.lb.an <- finite_rows(methods$data.lb.an)
  methods$data.lp.an <- finite_rows(methods$data.lp.an)

  alpha_true <- alpha_true_for(scenario$param, scenario$event, scenario$hypothesis)
  result <- do.call(cbind, lapply(data_names, function(name) {
    est_result(methods[[name]], alpha_true)
  }))
  colnames(result) <- method_order
  rownames(result) <- c("bias", "se", "acc", "coverage", "p")

  est <- do.call(cbind, lapply(data_names, function(name) methods[[name]]$estimate))
  se <- do.call(cbind, lapply(data_names, function(name) methods[[name]]$se))
  colnames(est) <- method_order
  colnames(se) <- method_order

  sid <- scenario_id(scenario)
  plot_estimates(est, alpha_true, scenario, file.path(plot_dir, paste0("ITC_est_", sid, ".png")))
  plot_se(se, scenario, file.path(plot_dir, paste0("ITC_se_", sid, ".png")))
  plot_metric(result, "coverage", 0.95, scenario, file.path(plot_dir, paste0("ITC_coverage_rate_", sid, ".png")))
  plot_metric(result, "acc", 1, scenario, file.path(plot_dir, paste0("ITC_accuracy_", sid, ".png")), lower = 0, upper = 2)

  long <- do.call(rbind, lapply(method_order, function(method) {
    data.frame(
      param = scenario$param,
      event = scenario$event,
      hypothesis = scenario$hypothesis,
      n = scenario$n,
      ess = scenario$ess,
      method = method,
      bias = result["bias", method],
      se = result["se", method],
      acc = result["acc", method],
      coverage = result["coverage", method],
      p = result["p", method],
      used_lbp = methods$used_lbp,
      used_rlp_file = methods$used_rlp_file,
      stringsAsFactors = FALSE
    )
  }))

  list(result = result, summary = long)
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

file_index <- parse_result_files(input_dir)
scenario_status <- find_complete_scenarios(file_index)
scenario_status <- filter_scenarios(scenario_status)
write.csv(
  scenario_status,
  file.path(plot_dir, "anchored_complete_scenarios_status.csv"),
  row.names = FALSE
)

complete_scenarios <- scenario_status[scenario_status$complete, ]
if (nrow(complete_scenarios) == 0) {
  stop("No complete scenarios found.")
}

summaries <- list()
for (i in seq_len(nrow(complete_scenarios))) {
  scenario <- complete_scenarios[i, ]
  message(sprintf(
    "[%d/%d] %s",
    i, nrow(complete_scenarios), scenario_id(scenario)
  ))
  summaries[[i]] <- analyze_scenario(scenario)$summary
}

summary_table <- do.call(rbind, summaries)
summary_table <- summary_table[order(
  summary_table$param, summary_table$event, summary_table$hypothesis,
  summary_table$n, summary_table$ess, summary_table$method
), ]

write.csv(
  summary_table,
  file.path(plot_dir, "anchored_results_summary.csv"),
  row.names = FALSE
)

message(sprintf("Complete scenarios analyzed: %d", nrow(complete_scenarios)))
message(sprintf("Plots and summary written to: %s", normalizePath(plot_dir, winslash = "/")))
