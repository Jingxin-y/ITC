## Independent RD generation only. Run from anchored_v1.
source("load_anchored.R", encoding = "UTF-8")
source_anchored_data()
design <- webster_rd_design(event = "common", ess_ratio = .5)
set.seed(1)
sim <- generate_anchored_data(n_AC = 1000, n_BC = 1000, design = design)
design$parameters
design$marginal_truth[, c("p_C", "p_A", "p_B", "RD_AC", "RD_BC", "RD_AB")]
design$rd_coefficient_truth # RD intercept/slope in fixed h_E units, NOT RR alpha
design$calibration
head(sim$data.IPD$data)
head(sim$data.AD.or$data)
## Use event="rare" to switch to calibrated rare-event risks.
## Existing comparison/configuration remains RR; this example does not fit RD models.
