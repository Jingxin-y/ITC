## Run from anchored_v1: Rscript data/example_webster_RR_scenarios.R
code_dir <- Sys.getenv("MAIC_CODE_DIR", unset = ".")
source(file.path(code_dir, "load_anchored.R"), encoding = "UTF-8")
source_anchored_data(code_dir)

## Full theoretical ESS/n=.15; only the EM group differs between AC and BC.
design <- webster_rr_scenario_design(
  event = "rare", ess_ratio = .15, ess_source = "em_driven", em_level = "high"
)
set.seed(1)
sim <- generate_anchored_data(n_AC = 200, n_BC = 200, design = design)
print(design$scenario)
print(design$relative_ess)  # full=.15, em=.15, non_em=1
print(design$marginal_truth)
print(sim$data.IPD$count)
print(sim$data.AD.or$count)

## Equivalent convenience entry, with seed handled inside the function:
same <- generate_webster_rr_scenario(
  n_AC = 200, n_BC = 200, event = "rare", ess_ratio = .15,
  ess_source = "em_driven", em_level = "high", seed = 1
)
stopifnot(identical(sim, same))

## Existing estimation entry: load estimation functions first.
## Choose each MAIC calibration set on the SAME generated dataset if comparing them.
## source_anchored(code_dir)
## scores <- make_brm_score_spec(paste0("z",1:3), paste0("z",4:6))
## result <- run.anchored(design, 200, 200, seed=1, score_spec=scores,
##   weight_covariates=paste0("z",1:6), model_weighting="unit",
##   bc_target="marginal", run_exact=TRUE)
