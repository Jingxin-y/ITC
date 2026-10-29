## From anchored_v1: source("RD/example_rd.R", encoding="UTF-8")
if (!exists("source_anchored_rd")) source("RD/load_anchored_rd.R", encoding = "UTF-8")
source_anchored_rd()
design <- webster_rd_scenario_design("common", .5, "diffuse", "high")
result <- run.anchored.rd(design, 500, seed = 19, run_exact = FALSE,
  brm_effect_model = "linear_atanh")
print(design$marginal_truth[, c("p_C", "p_A", "p_B", "RD_AC", "RD_BC", "RD_AB")])
print(result$estimates[, c("method", "estimate", "se", "lower", "upper", "status", "interval_status")])
