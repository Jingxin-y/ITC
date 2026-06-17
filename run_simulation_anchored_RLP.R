.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(doSNOW)
  library(doRNG)
  library(brm)
  library(MASS)
  library(sandwich)
})

param <- "RR"
hypothesis <- "null"
R <- 1000
first_r <- NA_integer_
last_r <- NA_integer_
result_dir <- Sys.getenv("result_dir", getwd())
code_dir <- Sys.getenv("code_dir", "/home/yanjin41/brmplus_simulation/compareITC/sandwich")
run_exact <- TRUE
include_gcomp <- TRUE
n_star <- 50
n_boot <- 30
skip_existing <- TRUE
ess_side <- "lower"
target_mean <- NA_real_
n_cov <- 3
ipd_prob <- c(0.8, 0.2, 0.8)
ad_prob <- NA_real_
direction <- c(-1, 1, -1)

argv <- commandArgs(TRUE)
if (length(argv) > 0) {
  for (i in seq_along(argv)) {
    eval(parse(text = argv[[i]]))
  }
}

if (is.na(first_r)) first_r <- 1L
if (is.na(last_r)) last_r <- R
seeds <- seq.int(first_r, last_r)

null_if_na <- function(x) {
  if (length(x) == 1 && is.na(x)) NULL else x
}

format_ess <- function(x) {
  sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

target_mean_arg <- null_if_na(target_mean)
ad_prob_arg <- null_if_na(ad_prob)
direction_arg <- null_if_na(direction)

scenarios <- do.call(rbind, lapply(c("common", "rare"), function(event_i) {
  expand.grid(
    event = event_i,
    n = c(100),
    ess_ratio = c(0.15, 0.25, 0.35, 0.5, 0.8),
    stringsAsFactors = FALSE
  )
}))

source_compareitc <- function(code_dir) {
  source(file.path(code_dir, "getProbScalarRR.R"))
  source(file.path(code_dir, "getProbScalarRD.R"))
  source(file.path(code_dir, "MyFunc.R"))
  source(file.path(code_dir, "1.1_MLE_Point.R"))
  source(file.path(code_dir, "1.2_MLE_Var.R"))
  source(file.path(code_dir, "1.2_MLE_Var_joint.R"))
  source(file.path(code_dir, "1_CallMLE.R"))
  source(file.path(code_dir, "CI_exact_fast_RLP.R"))
  source(file.path(code_dir, "generate_RR.R"))
  source(file.path(code_dir, "compare_RR_sanchored_RLP.R"))

  rcpp_exports <- file.path(dirname(dirname(code_dir)), "R", "RcppExports.R")
  if (file.exists(rcpp_exports)) {
    source(rcpp_exports)
  } else {
    warning("RcppExports.R was not found: ", rcpp_exports, call. = FALSE)
  }
}

source_compareitc(code_dir)

empty_anchored_result <- function() {
  out <- matrix(NA_real_, nrow = 5, ncol = 11)
  rownames(out) <- c("point.est", "se.est", "con.low", "con.up", "p.value")
  colnames(out) <- c(
    "brm", "brm_an", "LB_an", "LP_an", "RLP_an",
    "CMH_an", "CMH_sandwich_an", "brmad_an", "GC", "exact_an", "ad_exact_an"
  )
  out
}

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
ncores <- min(ncores, length(seeds), 101)
cl <- makeCluster(ncores, type = "SOCK")
on.exit(stopCluster(cl), add = TRUE)

clusterExport(cl, c("code_dir", "empty_anchored_result"), envir = environment())
clusterEvalQ(cl, {
  .libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))
  suppressPackageStartupMessages({
    library(brm)
    library(MASS)
    library(sandwich)
  })
  source(file.path(code_dir, "getProbScalarRR.R"))
  source(file.path(code_dir, "getProbScalarRD.R"))
  source(file.path(code_dir, "MyFunc.R"))
  source(file.path(code_dir, "1.1_MLE_Point.R"))
  source(file.path(code_dir, "1.2_MLE_Var.R"))
  source(file.path(code_dir, "1.2_MLE_Var_joint.R"))
  source(file.path(code_dir, "1_CallMLE.R"))
  source(file.path(code_dir, "CI_exact_fast_RLP.R"))
  source(file.path(code_dir, "generate_RR.R"))
  source(file.path(code_dir, "compare_RR_sanchored_RLP.R"))
  rcpp_exports <- file.path(dirname(dirname(code_dir)), "R", "RcppExports.R")
  if (file.exists(rcpp_exports)) source(rcpp_exports)
  NULL
})

registerDoSNOW(cl)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cat("Running all anchored null scenarios with joint RLP variance\n")
cat("param =", param, "hypothesis =", hypothesis, "R =", R,
    "seeds =", paste0(first_r, ":", last_r),
    "ncores =", ncores, "\n")
cat("run_exact =", run_exact, "include_gcomp =", include_gcomp,
    "n_star =", n_star, "n_boot =", n_boot,
    "skip_existing =", skip_existing, "\n")

run_one_scenario <- function(scenario) {
  event_i <- scenario$event
  n_i <- scenario$n
  ess_ratio_i <- scenario$ess_ratio
  ess_ratio_arg <- null_if_na(ess_ratio_i)
  ess_label <- format_ess(ess_ratio_i)

  out_file <- file.path(
    result_dir,
    paste0(
      "ITC_simulation_results_anchored_RLP_",
      param, "_", event_i, "_", hypothesis,
      "_n_", n_i, "_R_", R,
      "_ess_", ess_label,
      ".csv"
    )
  )

  if (isTRUE(skip_existing) && file.exists(out_file)) {
    cat("Skipping existing result:\n", out_file, "\n", sep = "")
    return(invisible(out_file))
  }

  cat("Running scenario: event =", event_i,
      "n =", n_i,
      "ess_ratio =", ess_label,
      "\n")

  result.mle <- foreach(
    r = seeds,
    .packages = c("brm", "MASS", "sandwich"),
    .export = c(
      "empty_anchored_result",
      "param", "hypothesis", "run_exact", "include_gcomp",
      "n_star", "n_boot", "ess_side", "target_mean_arg",
      "n_cov", "ipd_prob", "ad_prob_arg", "direction_arg",
      "n_i", "event_i", "ess_ratio_arg", "ess_label"
    )
  ) %dopar% {
    set.seed(r)
    error_message <- NA_character_

    r1 <- tryCatch(
      run.anchored(
        param,
        n_i,
        event_i,
        hypothesis,
        run_exact = run_exact,
        include_gcomp = include_gcomp,
        n_star = n_star,
        n_boot = n_boot,
        ess_ratio = ess_ratio_arg,
        ess_side = ess_side,
        target_mean = target_mean_arg,
        n_cov = n_cov,
        ipd_prob = ipd_prob,
        ad_prob = ad_prob_arg,
        direction = direction_arg
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        empty_anchored_result()
      }
    )

    list(
      r = r,
      estimate = r1[1, ],
      se = r1[2, ],
      low = r1[3, ],
      up = r1[4, ],
      p = r1[5, ],
      error = error_message
    )
  }

  result.all <- do.call(rbind, lapply(result.mle, function(x) {
    as.data.frame(x[c("r", "estimate", "se", "low", "up", "p")])
  }))
  write.csv(result.all, file = out_file, row.names = FALSE)

  errors <- do.call(rbind, lapply(result.mle, function(x) {
    if (is.na(x$error)) return(NULL)
    data.frame(
      event = event_i,
      n = n_i,
      ess_ratio = ess_label,
      r = x$r,
      error = x$error,
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(errors) && nrow(errors) > 0) {
    error_file <- sub("\\.csv$", "_errors.csv", out_file)
    write.csv(errors, file = error_file, row.names = FALSE)
    cat("Wrote seed-level errors:\n", error_file, "\n", sep = "")
    cat("First error: ", errors$error[1], "\n", sep = "")
  }

  cat("Wrote result:\n", out_file, "\n", sep = "")
  invisible(out_file)
}

for (i in seq_len(nrow(scenarios))) {
  cat("[", i, "/", nrow(scenarios), "] ", sep = "")
  run_one_scenario(scenarios[i, ])
}

cat("Finished all anchored RLP null scenarios.\n")
