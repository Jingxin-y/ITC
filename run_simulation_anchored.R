.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(doSNOW)
  library(doRNG)
  library(brm)
  library(MASS)
  library(sandwich)
})

param <- "RR"
event <- "common"
hypothesis <- "alternative"
n <- 50
R <- 400
first_r <- NA_integer_
last_r <- NA_integer_
result_dir <- Sys.getenv("result_dir", getwd())
code_dir <- Sys.getenv("code_dir", "/home/yanjin41/brmplus_simulation/compareITC/sandwich")
run_exact <- TRUE
include_gcomp <- TRUE
n_star <- 500
n_boot <- 300
ess_ratio <- NA_real_
ess_side <- "lower"
target_mean <- NA_real_
n_cov <- 3
ipd_prob <- c(0.7, 0.3, 0.7)
ad_prob <- NA_real_
direction <- c(-1, 1, -1)

argv <- commandArgs(TRUE)
if (length(argv) > 0) {
  for (i in seq_along(argv)) {
    eval(parse(text = argv[[i]]))
  }
}

if (is.na(first_r)) first_r <- R - 199
if (is.na(last_r)) last_r <- R
seeds <- seq.int(first_r, last_r)
null_if_na <- function(x) {
  if (length(x) == 1 && is.na(x)) NULL else x
}
ess_ratio_arg <- null_if_na(ess_ratio)
target_mean_arg <- null_if_na(target_mean)
ad_prob_arg <- null_if_na(ad_prob)
direction_arg <- null_if_na(direction)

source_compareitc <- function(code_dir) {
  source(file.path(code_dir, "getProbScalarRR.R"))
  source(file.path(code_dir, "getProbScalarRD.R"))
  source(file.path(code_dir, "MyFunc.R"))
  source(file.path(code_dir, "1.1_MLE_Point.R"))
  source(file.path(code_dir, "1.2_MLE_Var.R"))
  source(file.path(code_dir, "1.2_MLE_Var_joint.R"))
  source(file.path(code_dir, "1_CallMLE.R"))
  source(file.path(code_dir, "CI_exact_fast.R"))
  source(file.path(code_dir, "generate_RR.R"))
  source(file.path(code_dir, "compare_RR_anchored.R"))

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

clusterExport(cl, "code_dir", envir = environment())
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
  source(file.path(code_dir, "CI_exact_fast.R"))
  source(file.path(code_dir, "generate_RR.R"))
  source(file.path(code_dir, "compare_RR_anchored.R"))
  rcpp_exports <- file.path(dirname(dirname(code_dir)), "R", "RcppExports.R")
  if (file.exists(rcpp_exports)) source(rcpp_exports)
  NULL
})

registerDoSNOW(cl)

cat("Running anchored-only simulation\n")
cat("param =", param, "n =", n, "R =", R,
    "event =", event, "hypothesis =", hypothesis, "\n")
cat("seeds =", paste(seeds, collapse = ","), "\n")
cat("run_exact =", run_exact, "include_gcomp =", include_gcomp,
    "n_star =", n_star, "n_boot =", n_boot, "\n")
cat("ESS controls: ess_ratio =", ifelse(is.null(ess_ratio_arg), "NULL", ess_ratio_arg),
    "ess_side =", ess_side,
    "target_mean =", ifelse(is.null(target_mean_arg), "NULL", target_mean_arg),
    "n_cov =", n_cov,
    "ipd_prob =", paste(ipd_prob, collapse = ","),
    "ad_prob =", ifelse(is.null(ad_prob_arg), "NULL", paste(ad_prob_arg, collapse = ",")),
    "direction =", ifelse(is.null(direction_arg), "NULL", paste(direction_arg, collapse = ",")),
    "\n")

result.mle <- foreach(
  r = seeds,
  .packages = c("brm", "MASS", "sandwich")
) %dopar% {
  set.seed(r)

  r1 <- tryCatch(
    run.anchored(
      param,
      n,
      event,
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
      warning("seed ", r, " failed: ", conditionMessage(e), call. = FALSE)
      empty_anchored_result()
    }
  )

  list(
    r = r,
    estimate = r1[1, ],
    se = r1[2, ],
    low = r1[3, ],
    up = r1[4, ],
    p = r1[5, ]
  )
}

result.all <- do.call(rbind, lapply(result.mle, as.data.frame))


dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(
  result_dir,
  paste0(
    "ITC_simulation_results_anchored_",
    param, "_", event, "_", hypothesis,
    "_n_", n, "_R_", R,
    "_ess_", ess_ratio,
    ".csv"
  )
)

write.csv(result.all, file = out_file, row.names = FALSE)
cat("Wrote anchored-only results to:\n", out_file, "\n", sep = "")
