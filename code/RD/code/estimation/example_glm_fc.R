## Run from anchored_v1 after loading the estimation code.
source("load_anchored.R", encoding = "UTF-8")
source_anchored()

example_glm_fc <- function(event = "rare", model_weighting = "unit", n = 500L, seed = 811L) {
  result <- run.anchored(webster_rr_design(event = event, ess_ratio = .5),
    n_AC = n, n_BC = n, seed = seed, score_spec = make_brm_score_spec(),
    weight_covariates = paste0("z", 1:6), model_weighting = model_weighting,
    run_exact = FALSE, return_full = TRUE, return_data = TRUE)
  print(subset(result$estimates,
    method %in% c("LB", "LP", "RLP", "RLB", anchored_glm_fc_method_names()) &
      term == "x", select = c(method, estimate, se, lower, upper, status)))
  invisible(result)
}

## result <- example_glm_fc()
## result <- example_glm_fc(model_weighting = "maic")
## Coefficient estimates are on the log-RR scale. A BC-target marginal estimate
## additionally standardizes predictions using result$weights$weight.
## source("estimation/cluster_all_methods.R", encoding = "UTF-8")
## marginal <- maic_wald_rows(result)
## subset(marginal, method %in% anchored_glm_fc_method_names())
