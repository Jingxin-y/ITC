## 旧 DGM 的批量比较入口；尚未切换至联合 RR-GOP / effect modifier 分析。
## 在 anchored_v1 工作目录运行；BC 方法参考 brmplus/development。

suppressPackageStartupMessages({
  library(doSNOW)
  library(doRNG)
  library(MASS)
})

param <- "RR"
event <- "common"
hypothesis <- "alternative"
n <- 50
R <- 400
first_r <- NA_integer_
last_r <- NA_integer_
result_dir <- Sys.getenv("result_dir", "results")
code_dir <- Sys.getenv("code_dir", ".")
run_exact <- TRUE
exact_seed_offset <- 1000000L
bc_control <- list()
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
  source(file.path(code_dir, "load_anchored.R"))
  source_anchored(code_dir, legacy = TRUE)
}

source_compareitc(code_dir)

empty_anchored_result <- function() {
  out <- matrix(NA_real_, nrow = 5, ncol = 13)
  rownames(out) <- c("point.est", "se.est", "con.low", "con.up", "p.value")
  colnames(out) <- c(
    "brm", "brm_an", "LB_an", "LP_an", "RLP_an",
    "CMH_an", "CMH_sandwich_an", "brmad_an", "brm-BC", "brm-FC-BC", "brm_b-BC", "RLB_an", "RLB_sandwich_an"
  )
  out
}

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
ncores <- min(ncores, length(seeds), 101)
cl <- makeCluster(ncores, type = "SOCK")

clusterExport(cl, "code_dir", envir = environment())
clusterEvalQ(cl, {
  source(file.path(code_dir, "load_anchored.R"))
  source_anchored(code_dir, legacy = TRUE)
  NULL
})

registerDoSNOW(cl)

cat("Running anchored-only simulation\n")
cat("param =", param, "n =", n, "R =", R,
    "event =", event, "hypothesis =", hypothesis, "\n")
cat("seeds =", paste(seeds, collapse = ","), "\n")
cat("run_exact =", run_exact, "exact_seed_offset =", exact_seed_offset, "\n")
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
  .packages = c("MASS")
) %dopar% {
  set.seed(r)

  r1 <- tryCatch(
    run.anchored(
      param,
      n,
      event,
      hypothesis,
      run_exact = run_exact,
      ess_ratio = ess_ratio_arg,
      ess_side = ess_side,
      target_mean = target_mean_arg,
      n_cov = n_cov,
      ipd_prob = ipd_prob,
      ad_prob = ad_prob_arg,
      direction = direction_arg,
      seed = r, exact_seed_offset = exact_seed_offset, bc_control = bc_control
    ),
    error = function(e) {
      warning("seed ", r, " failed: ", conditionMessage(e), call. = FALSE)
      empty_anchored_result()
    }
  )

  list(
    r = r,
    method = colnames(r1),
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

stopCluster(cl)
