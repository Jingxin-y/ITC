## Run from anchored_v1: Rscript RD/estimation/run_simulation_rd.R RD/simulation_config_rd.R RD/results/local
design <- webster_rd_scenario_design(event = "common", ess_ratio = .5,
  ess_source = "diffuse", em_level = "high")
simulation_args <- list(design = design, n_AC = 1000, n_BC = 1000, seeds = 1:10,
  run_exact = FALSE, weight_covariates = paste0("z", 1:6), model_weighting = "unit",
  score_standardization = "bc_estimated", brm_effect_model = "linear_atanh",
  fit_control = list(), bc_control = list())
