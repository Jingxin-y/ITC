## From anchored_v1: Rscript estimation/tests/check_webster_rr_scenario_pipeline.R
code_dir <- Sys.getenv("MAIC_CODE_DIR", unset = ".")
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
source_anchored(code_dir)
design <- webster_rr_scenario_design("common", .5, "nonem_driven", "high")
scores <- make_brm_score_spec(paste0("z", 1:3), paste0("z", 4:6))
result <- suppressWarnings(run.anchored(
  design, n_AC = 1000, n_BC = 1000, seed = 811,
  score_spec = scores, weight_covariates = paste0("z", 1:6),
  model_weighting = "unit", bc_target = "marginal", run_exact = FALSE, return_data = TRUE
))
stopifnot(identical(result$simulation$design, design),
          identical(result$variables$model_weighting, "unit"),
          identical(result$variables$standardization_weighting, "maic"),
          nrow(result$estimates) > 0L, nrow(result$reference_truth) > 0L,
          all(is.finite(result$reference_truth$truth)),
          nrow(result$data$data.IPD$data) == 1000L)
lp <- subset(result$estimates, method == "LP" & component == "treatment")
stopifnot(nrow(lp) == 2L, all(lp$status == "available"), all(is.finite(lp$estimate)))
cat("PASS: new DGP integrates with existing estimates, unit model weights, MAIC targets and replicate coefficient truth\n")
