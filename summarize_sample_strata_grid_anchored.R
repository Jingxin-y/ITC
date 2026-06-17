## Batch reconstruction of anchored simulation samples.
##
## Defaults scan:
##   n = 50,100,200,500
##   event = common,rare
##   hypothesis = alternative,null
##   ess_ratio = 0.15,0.25,0.35,0.5,0.8
##   R = 1000, seeds = 1:R
##
## Outputs:
##   *_by_seed.csv: one row per seed x group x covariate stratum
##   *_summary.csv: across-seed summary for that scenario
##   ITC_sample_strata_anchored_grid_summary.csv: combined summary

param <- "RR"
R <- 1000
first_r <- 1L
last_r <- R
n_values <- c(50, 100, 200, 500)
events <- c("common", "rare")
hypotheses <- c("alternative", "null")
ess_ratios <- c(0.15, 0.25, 0.35, 0.5, 0.8)
result_dir <- Sys.getenv("result_dir", getwd())
output_dir <- Sys.getenv("output_dir", file.path(result_dir, "sample_strata_grid"))
code_dir <- Sys.getenv("code_dir", getwd())
ess_side <- "lower"
target_mean <- NA_real_
n_cov <- 3
ipd_prob <- c(0.8, 0.2, 0.8)
ad_prob <- NA_real_
direction <- c(-1, 1, -1)
include_empty_strata <- TRUE
write_by_seed <- TRUE
write_scenario_summary <- TRUE

argv <- commandArgs(TRUE)
if (length(argv) > 0) {
  for (i in seq_along(argv)) {
    eval(parse(text = argv[[i]]))
  }
}

## Accept the common typo/variant from notes or command lines.
if (exists("rss_ratio", inherits = FALSE)) {
  ess_ratios <- get("rss_ratio", inherits = FALSE)
}
if (exists("rss_ratios", inherits = FALSE)) {
  ess_ratios <- get("rss_ratios", inherits = FALSE)
}
if (exists("ess_ratio", inherits = FALSE)) {
  ess_ratios <- get("ess_ratio", inherits = FALSE)
}

null_if_na <- function(x) {
  if (length(x) == 1 && is.na(x)) NULL else x
}

format_value <- function(x) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) {
    return("NA")
  }
  sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

source_required <- function(code_dir) {
  source(file.path(code_dir, "MyFunc.R"))
  source(file.path(code_dir, "getProbScalarRR.R"))
  source(file.path(code_dir, "generate_RR.R"))
  source(file.path(code_dir, "compare_RR_anchored.R"))
}

covariate_key <- function(data, z_cols) {
  if (length(z_cols) == 0) {
    return(rep("all", nrow(data)))
  }
  do.call(paste, c(data[z_cols], sep = ":"))
}

make_strata_frame <- function(data, z_cols, include_empty_strata) {
  if (length(z_cols) == 0) {
    return(data.frame(stratum = "all", stringsAsFactors = FALSE))
  }

  if (isTRUE(include_empty_strata)) {
    out <- expand.grid(
      rep(list(0:1), length(z_cols)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(out) <- z_cols
  } else {
    out <- unique(data[z_cols])
    out <- out[do.call(order, out), , drop = FALSE]
  }
  out$stratum <- covariate_key(out, z_cols)
  out
}

weighted_ess <- function(w) {
  if (length(w) == 0 || sum(w) <= 0) return(0)
  sum(w)^2 / sum(w^2)
}

summarize_sample <- function(sample_data, seed, group, scenario,
                             include_empty_strata, weights = NULL) {
  data <- sample_data$data
  z_cols <- covariate_names(data)
  strata_frame <- make_strata_frame(data, z_cols, include_empty_strata)
  observed_key <- covariate_key(data, z_cols)

  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(data)) {
    stop("Length of weights must equal nrow(data).")
  }

  rows <- lapply(seq_len(nrow(strata_frame)), function(i) {
    key_i <- strata_frame$stratum[i]
    idx <- observed_key == key_i
    x0 <- idx & data$x == 0
    x1 <- idx & data$x == 1

    w <- weights[idx]
    w0 <- weights[x0]
    w1 <- weights[x1]

    base <- data.frame(
      scenario_n = scenario$n,
      event = scenario$event,
      hypothesis = scenario$hypothesis,
      ess_ratio = scenario$ess_ratio,
      r = seed,
      group = group,
      stringsAsFactors = FALSE
    )
    if (length(z_cols) > 0) {
      base <- cbind(base, strata_frame[i, z_cols, drop = FALSE])
    }

    weighted_n <- sum(w)
    weighted_y1 <- sum(weights[idx] * (data$y[idx] == 1))
    weighted_y0 <- sum(weights[idx] * (data$y[idx] == 0))

    cbind(
      base,
      data.frame(
        stratum = key_i,
        stratum_n = sum(idx),
        y1 = sum(data$y[idx] == 1),
        y0 = sum(data$y[idx] == 0),
        x0_n = sum(x0),
        x0_y1 = sum(data$y[x0] == 1),
        x0_y0 = sum(data$y[x0] == 0),
        x1_n = sum(x1),
        x1_y1 = sum(data$y[x1] == 1),
        x1_y0 = sum(data$y[x1] == 0),
        weighted_n = weighted_n,
        weighted_y1 = weighted_y1,
        weighted_y0 = weighted_y0,
        weighted_event_rate = if (weighted_n > 0) weighted_y1 / weighted_n else NA_real_,
        weighted_ess = weighted_ess(w),
        x0_weighted_n = sum(w0),
        x0_weighted_y1 = sum(weights[x0] * (data$y[x0] == 1)),
        x0_weighted_y0 = sum(weights[x0] * (data$y[x0] == 0)),
        x1_weighted_n = sum(w1),
        x1_weighted_y1 = sum(weights[x1] * (data$y[x1] == 1)),
        x1_weighted_y0 = sum(weights[x1] * (data$y[x1] == 0)),
        stringsAsFactors = FALSE
      )
    )
  })

  do.call(rbind, rows)
}

summarize_across_seeds <- function(by_seed) {
  split_cols <- c(
    "scenario_n", "event", "hypothesis", "ess_ratio",
    "group", covariate_names(by_seed), "stratum"
  )
  pieces <- split(
    by_seed,
    interaction(by_seed[split_cols], drop = TRUE, lex.order = TRUE),
    drop = TRUE
  )

  out <- lapply(pieces, function(df) {
    id <- df[1, split_cols, drop = FALSE]
    total_n <- sum(df$stratum_n)
    total_weighted_n <- sum(df$weighted_n)

    data.frame(
      id,
      replications = nrow(df),
      total_n = total_n,
      total_y1 = sum(df$y1),
      total_y0 = sum(df$y0),
      event_rate = if (total_n > 0) sum(df$y1) / total_n else NA_real_,
      total_weighted_n = total_weighted_n,
      total_weighted_y1 = sum(df$weighted_y1),
      total_weighted_y0 = sum(df$weighted_y0),
      weighted_event_rate = if (total_weighted_n > 0) {
        sum(df$weighted_y1) / total_weighted_n
      } else {
        NA_real_
      },
      mean_n = mean(df$stratum_n),
      sd_n = stats::sd(df$stratum_n),
      min_n = min(df$stratum_n),
      max_n = max(df$stratum_n),
      mean_y1 = mean(df$y1),
      mean_y0 = mean(df$y0),
      mean_weighted_n = mean(df$weighted_n),
      sd_weighted_n = stats::sd(df$weighted_n),
      min_weighted_n = min(df$weighted_n),
      max_weighted_n = max(df$weighted_n),
      mean_weighted_ess = mean(df$weighted_ess),
      empty_replications = sum(df$stratum_n == 0),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[do.call(order, out[split_cols]), , drop = FALSE]
}

scenario_id <- function(scenario) {
  paste0(
    "RR_", scenario$event, "_", scenario$hypothesis,
    "_n_", scenario$n,
    "_R_", R,
    "_ess_", format_value(scenario$ess_ratio)
  )
}

run_scenario <- function(scenario, seeds) {
  ess_ratio_arg <- null_if_na(scenario$ess_ratio)
  target_mean_arg <- null_if_na(target_mean)
  ad_prob_arg <- null_if_na(ad_prob)
  direction_arg <- null_if_na(direction)

  by_seed <- lapply(seeds, function(seed) {
    set.seed(seed)

    data_ipd <- data.generation(
      param, "IPD", scenario$n, scenario$event, scenario$hypothesis,
      n_cov = n_cov,
      ipd_prob = ipd_prob
    )
    data_ad_or <- data.generation(
      param, "AD", scenario$n, scenario$event, scenario$hypothesis,
      ess_ratio = ess_ratio_arg,
      ess_side = ess_side,
      target_mean = target_mean_arg,
      n_cov = n_cov,
      ipd_prob = ipd_prob,
      ad_prob = ad_prob_arg,
      direction = direction_arg
    )

    data_ad <- get.aggre(data_ad_or)
    anchored_weight <- get.weight(data_ipd, data_ad$anchored.mean.x)$weight

    do.call(rbind, list(
      summarize_sample(
        data_ipd, seed, "IPD", scenario, include_empty_strata
      ),
      summarize_sample(
        data_ad_or, seed, "AD", scenario, include_empty_strata
      ),
      summarize_sample(
        data_ipd, seed, "IPD_REWEIGHTED", scenario,
        include_empty_strata, weights = anchored_weight
      )
    ))
  })

  do.call(rbind, by_seed)
}

if (!dir.exists(output_dir)) {
  ok <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(ok) || !dir.exists(output_dir)) {
    stop("Cannot create output_dir: ", output_dir)
  }
}

source_required(code_dir)

seeds <- seq.int(as.integer(first_r), as.integer(last_r))
scenarios <- expand.grid(
  n = n_values,
  event = events,
  hypothesis = hypotheses,
  ess_ratio = ess_ratios,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

cat("Summarizing anchored sample strata grid\n")
cat("scenarios =", nrow(scenarios), "\n")
cat("R label =", R, "seeds =", paste(range(seeds), collapse = ":"), "\n")
cat("output_dir =", output_dir, "\n")

manifest <- data.frame()
summaries <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, ]
  sid <- scenario_id(scenario)
  cat(sprintf("[%d/%d] %s\n", i, nrow(scenarios), sid))

  by_seed <- run_scenario(scenario, seeds)
  summary <- summarize_across_seeds(by_seed)
  summaries[[i]] <- summary

  by_seed_file <- file.path(output_dir, paste0("ITC_sample_strata_anchored_", sid, "_by_seed.csv"))
  summary_file <- file.path(output_dir, paste0("ITC_sample_strata_anchored_", sid, "_summary.csv"))

  if (isTRUE(write_by_seed)) {
    write.csv(by_seed, by_seed_file, row.names = FALSE)
  } else {
    by_seed_file <- NA_character_
  }
  if (isTRUE(write_scenario_summary)) {
    write.csv(summary, summary_file, row.names = FALSE)
  } else {
    summary_file <- NA_character_
  }

  manifest <- rbind(
    manifest,
    data.frame(
      scenario_n = scenario$n,
      event = scenario$event,
      hypothesis = scenario$hypothesis,
      ess_ratio = scenario$ess_ratio,
      R = R,
      first_r = min(seeds),
      last_r = max(seeds),
      by_seed_file = by_seed_file,
      summary_file = summary_file,
      stringsAsFactors = FALSE
    )
  )
}

combined_summary <- do.call(rbind, summaries)
combined_summary <- combined_summary[do.call(order, combined_summary[
  c("scenario_n", "event", "hypothesis", "ess_ratio", "group", covariate_names(combined_summary), "stratum")
]), , drop = FALSE]

combined_summary_file <- file.path(
  output_dir,
  "ITC_sample_strata_anchored_grid_summary.csv"
)
manifest_file <- file.path(
  output_dir,
  "ITC_sample_strata_anchored_grid_manifest.csv"
)

write.csv(combined_summary, combined_summary_file, row.names = FALSE)
write.csv(manifest, manifest_file, row.names = FALSE)

cat("Wrote combined summary to:\n", combined_summary_file, "\n", sep = "")
cat("Wrote manifest to:\n", manifest_file, "\n", sep = "")
