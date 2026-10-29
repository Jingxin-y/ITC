## Direct additive risks. OP/GOP inversion is never used to generate data.
webster_rd_scenario_targets <- function(event = "common") {
  event <- webster_event(event)
  target <- webster_event_targets(event)["BC", ]
  target["p_A"] <- target["p_B"] + if (event == "common") .14 else .05
  target
}

webster_rd_scenario_parameters <- function(event = "common", em_level = c("low", "high")) {
  event <- webster_event(event); em_level <- match.arg(em_level)
  target <- webster_rd_scenario_targets(event)
  preset <- if (event == "common") c(beta_E = -0.1684896018658,
    beta_P = 0.0164089587479, lambda_high = 0.1809727741321) else
    c(beta_E = -0.01999854953322, beta_P = 0.00203773239874,
      lambda_high = 0.02147941217947)
  lambda <- unname(preset["lambda_high"]) / if (em_level == "low") 4 else 1
  list(baseline_risk = unname(target["p_C"] - .9 * sum(preset[c("beta_E", "beta_P")])),
    beta_E = unname(preset["beta_E"]), beta_P = unname(preset["beta_P"]),
    rd_AC = unname(target["p_A"] - target["p_C"] - .9 * lambda / sqrt(.63)),
    rd_BC = unname(target["p_B"] - target["p_C"] - .45 * lambda / sqrt(.63)),
    em_rd_AC = lambda / sqrt(.63), em_rd_BC = .5 * lambda / sqrt(.63))
}

webster_rd_scenario_design <- function(event = "common", ess_ratio = .5,
    ess_source = c("diffuse", "em_driven", "nonem_driven"),
    em_level = c("low", "high"), parameters = list()) {
  event <- webster_event(event); ess_source <- match.arg(ess_source); em_level <- match.arg(em_level)
  population <- webster_rr_scenario_population(ess_ratio, ess_source)
  p <- webster_rd_scenario_parameters(event, em_level)
  if (!is.list(parameters) || (length(parameters) && (is.null(names(parameters)) ||
      anyDuplicated(names(parameters)) || any(!names(parameters) %in% names(p)))))
    stop("parameters must be a named list of additive-risk coefficients.")
  p[names(parameters)] <- parameters
  design <- do.call(webster_rd_design, c(list(event = event, q = population$q), p))
  target <- webster_rd_scenario_targets(event); achieved <- design$marginal_truth["BC", names(target)]
  error <- max(abs(achieved - target))
  if (!length(parameters) && error > 1e-11) stop("RD BC risk calibration failed.")
  margin <- if (event == "common") .05 else .005
  risks <- as.matrix(design$risk_table[c("p_C", "p_A", "p_B")])
  if (!length(parameters) && any(risks < margin - 1e-10 | risks > 1 - margin + 1e-10))
    stop("RD preset is outside its calibrated probability envelope.")
  design$dgm <- "webster2024_rd_ess_em_v2"
  design$ess_source <- ess_source; design$em_level <- em_level
  design$relative_ess <- population$relative_ess
  design$scenario <- list(version = 2L, ess_source = ess_source, em_level = em_level,
    nominal_full_ess_ratio = ess_ratio, log_ess_share_em = population$log_ess_share_em,
    lambda_AC = sqrt(.63) * p$em_rd_AC, lambda_BC = sqrt(.63) * p$em_rd_BC,
    rho = .5, reference_score_sd = sqrt(.63),
    target_RD_AB_BC = if (event == "common") .14 else .05,
    custom_parameters = parameters)
  design$calibration <- list(target_population = "BC", target_risks = target,
    achieved_risks = achieved, max_abs_error = error, matches_target_BC = error < 1e-11,
    reference_RR_BC_risks = webster_event_targets(event)["BC", ],
    intended_A_risk_shift = unname(target["p_A"] - webster_event_targets(event)["BC", "p_A"]),
    probability_envelope = c(margin, 1 - margin), overridden_parameters = names(parameters),
    reference_ess = c(.15, .25, .5, .8, 1),
    scope = "RD v2: BC A-B RD is .14 common / .05 rare; C/B risks retain RR targets; A risk is shifted uniformly from RD v1.")
  design$source$adaptation <- "Direct additive-risk Webster adaptation; three EMs and three prognostic non-EMs."
  design
}

generate_webster_rd_scenario <- function(n_AC, n_BC = n_AC, event = "common", ess_ratio = .5,
    ess_source = "diffuse", em_level = "high", seed = NULL, parameters = list()) {
  rd_validate_sizes(n_AC, n_BC)
  if (!is.null(seed)) { rd_validate_seeds(seed); set.seed(seed) }
  generate_anchored_data(n_AC, n_BC,
    webster_rd_scenario_design(event, ess_ratio, ess_source, em_level, parameters))
}

rd_calibration_table <- function() {
  rows <- list()
  for (event in c("common", "rare")) for (level in c("low", "high"))
    for (src in c("diffuse", "em_driven", "nonem_driven")) for (ess in c(.15, .25, .5, .8, 1)) {
      rd <- webster_rd_scenario_design(event, ess, src, level)
      rr <- webster_rr_scenario_design(event, ess, src, level)
      for (pop in c("AC", "BC")) for (arm in c("p_C", "p_A", "p_B"))
        rows[[length(rows) + 1L]] <- data.frame(event = event, em_level = level,
          ess_source = src, ess_ratio = ess, population = pop, arm = arm,
          preset_version = 2L,
          RD_risk = rd$marginal_truth[pop, arm], RR_risk = rr$marginal_truth[pop, arm],
          intended_risk_shift = if (arm == "p_A") rd$calibration$intended_A_risk_shift else 0,
          difference = rd$marginal_truth[pop, arm] - rr$marginal_truth[pop, arm])
    }
  do.call(rbind, rows)
}

rd_population_truth_table <- function(ess_ratios = c(.15, .25, .5, .8, 1)) {
  rows <- list()
  for (event in c("common", "rare")) for (level in c("low", "high"))
    for (src in c("diffuse", "em_driven", "nonem_driven")) for (ess in ess_ratios) {
      d <- webster_rd_scenario_design(event, ess, src, level)
      for (pop in c("AC", "BC")) rows[[length(rows) + 1L]] <- data.frame(
        event = event, em_level = level, ess_source = src, ess_ratio = ess,
        population = pop, preset_version = 2L,
        as.list(d$marginal_truth[pop, c("p_A", "p_B", "p_C", "RD_AC", "RD_BC", "RD_AB")]))
    }
  do.call(rbind, rows)
}
