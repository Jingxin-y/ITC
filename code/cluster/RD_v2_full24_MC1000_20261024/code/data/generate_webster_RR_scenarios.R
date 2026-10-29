## ESS-source / EM-strength extension of this project's direct-risk RR DGM.
## Load with source_anchored_data() or source_anchored().
## This is a project-specific design, not a reproduction of the paper's grid.
## BC means are fixed; AC means vary with the source population.

webster_rr_scenario_population <- function(
    ess_ratio = .5, ess_source = c("diffuse", "em_driven", "nonem_driven")) {
  ess_source <- match.arg(ess_source)
  if (!is.numeric(ess_ratio) || length(ess_ratio) != 1L ||
      !is.finite(ess_ratio) || ess_ratio <= 0 || ess_ratio > 1)
    stop("ess_ratio must be a finite scalar in (0,1].")
  ess_ratio <- unname(ess_ratio)

  e <- c(.7, .3, .7, .3, .7, .3)
  share <- switch(ess_source, diffuse = .5, em_driven = 1, nonem_driven = 0)
  group_targets <- c(em = ess_ratio^share, non_em = ess_ratio^(1 - share))
  q_E <- ipd_prob_for_ess_theory(group_targets["em"], e, c(-1, 1, -1, 0, 0, 0))
  q_P <- ipd_prob_for_ess_theory(group_targets["non_em"], e, c(0, 0, 0, 1, -1, 1))
  population <- webster_population(q = c(q_E[1:3], q_P[4:6]))
  relative_ess <- c(full = population$ess_ratio,
                   em = rel_ess_binary(population$q[1:3], population$e[1:3]),
                   non_em = rel_ess_binary(population$q[4:6], population$e[4:6]))
  target <- c(full = ess_ratio, group_targets)
  if (max(abs(relative_ess - target)) > 1e-9)
    stop("ESS calibration did not reach the requested population targets.")
  c(population, list(ess_source = ess_source, nominal_full_ess_ratio = ess_ratio,
                    log_ess_share_em = share, relative_ess = relative_ess))
}

webster_rr_scenario_parameters <- function(event = "common", em_level = c("low", "high")) {
  event <- webster_event(event)
  em_level <- match.arg(em_level)
  lambda <- switch(em_level, low = .15, high = .60)
  rho <- .5
  sigma <- sqrt(.63)  # Fixed BC population SD, never a replicate's estimated SD.
  g_E <- .55
  g_P <- 1.08
  target <- webster_event_targets(event)["BC", ]
  moment <- function(t) (.7 + .3 * t)^3
  em_AC <- exp(lambda / sigma)
  em_BC <- exp(rho * lambda / sigma)

  list(
    baseline_risk = unname(target["p_C"] / (moment(g_E) * moment(g_P))),
    rr_AC = unname((target["p_A"] / target["p_C"]) * moment(g_E) / moment(g_E * em_AC)),
    rr_BC = unname((target["p_B"] / target["p_C"]) * moment(g_E) / moment(g_E * em_BC)),
    em_rr_AC = em_AC, em_rr_BC = em_BC,
    baseline_rr_E = g_E, baseline_rr_P = g_P
  )
}

webster_rr_scenario_design <- function(
    event = "common", ess_ratio = .5,
    ess_source = c("diffuse", "em_driven", "nonem_driven"),
    em_level = c("low", "high")) {
  event <- webster_event(event)
  ess_source <- match.arg(ess_source)
  em_level <- match.arg(em_level)
  population <- webster_rr_scenario_population(ess_ratio, ess_source)
  parameters <- webster_rr_scenario_parameters(event, em_level)
  design <- do.call(webster_rr_design,
                    c(list(event = event, q = population$q), parameters))

  target <- webster_event_targets(event)["BC", c("p_C", "p_A", "p_B")]
  achieved <- design$marginal_truth["BC", names(target)]
  error <- max(abs(achieved - target))
  if (error > 1e-11) stop("Risk calibration did not preserve the BC population means.")
  risk <- as.matrix(design$risk_table[c("p_C", "p_A", "p_B")])
  if (any(!is.finite(risk)) || any(risk <= 0 | risk >= 1))
    stop("Scenario risks must be strictly between zero and one; no clipping is used.")

  lambda_AC <- unname(design$coefficient_truth["AC", "lambda"])
  lambda_BC <- unname(design$coefficient_truth["BC", "lambda"])
  design$dgm <- "webster2024_rr_ess_em_v1"
  design$ess_source <- ess_source
  design$em_level <- em_level
  design$relative_ess <- population$relative_ess
  design$scenario <- list(
    version = 1L, ess_source = ess_source, em_level = em_level,
    nominal_full_ess_ratio = ess_ratio, log_ess_share_em = population$log_ess_share_em,
    lambda_AC = lambda_AC, lambda_BC = lambda_BC, lambda_AB = lambda_AC - lambda_BC,
    rho = .5, reference_score_sd = sqrt(.63)
  )
  ## Replace the old calibration metadata: this extension only matches BC.
  design$calibration <- list(
    target_population = "BC", target_covariate_prob = design$e,
    target_risks = target, achieved_risks = achieved,
    max_abs_error = error, matches_target_BC = TRUE,
    calibrated_by = c("event", "em_level"),
    scope = "BC population C/A/B means fixed analytically; AC means unconstrained. Parameters do not depend on ESS source, ESS, sample size or seed."
  )
  design$source$adaptation <- paste(
    "Project extension: three EMs and three prognostic non-EMs; independent Bernoulli covariates;",
    "pure EM/non-EM or diffuse population shifts; low/high EM; fixed BC means; balanced allocation."
  )
  design
}

## n_AC and n_BC are TOTAL trial sizes; each observed arm contains n/2 people.
## Returned data.IPD / data.AD.or / data.AD / design match generate_anchored_data().
generate_webster_rr_scenario <- function(
    n_AC, n_BC = n_AC, event = "common", ess_ratio = .5,
    ess_source = c("diffuse", "em_driven", "nonem_driven"),
    em_level = c("low", "high"), seed = NULL) {
  for (size in list(n_AC, n_BC)) {
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
        size < 2 || size %% 2 != 0 || size > .Machine$integer.max)
      stop("n_AC and n_BC must be positive even integers within R's integer range.")
  }
  design <- webster_rr_scenario_design(event, ess_ratio, match.arg(ess_source), match.arg(em_level))
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
      stop("seed must be a nonnegative integer within R's integer range.")
    set.seed(seed)
  }
  generate_anchored_data(n_AC, n_BC, design)
}
