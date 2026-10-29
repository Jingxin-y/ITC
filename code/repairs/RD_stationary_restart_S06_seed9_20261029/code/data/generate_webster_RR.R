## Direct log-risk DGM inspired by the RR sensitivity analysis of
## Webster-Clark et al. (2024), DOI: 10.1093/aje/kwae048.
## Author code at revision 6d495d97b5002d16e5ae50bfd2f33c1a58f3dddc:
## https://github.com/mawcpharmdphd/ext_var_sel/blob/6d495d97b5002d16e5ae50bfd2f33c1a58f3dddc/Sensitivity_analyses/simulating_variable_selection_RR.sas
## This is an ADAPTATION: retain this project's THREE EMs (z1:z3), THREE
## non-EMs (z4:z6), score directions, BC distribution and ESS calibration.
## The paper instead uses two EMs, two prognostic non-EMs and two null variables.
## The B/C arm and rare-event preset are also project extensions.

## All coefficient inputs below are positive risk MULTIPLIERS on RAW scores.
## p_C = baseline_risk * baseline_rr_E^S_E * baseline_rr_P^S_P
## p_A = p_C * rr_AC * em_rr_AC^S_E
## p_B = p_C * rr_BC * em_rr_BC^S_E
## No OP/GOP outcome model, outcome root solving, or probability clipping is used.
## Defaults match the old common/rare marginal risks in BOTH populations
## at ESS=.5. The calibration is fixed when the population/ESS is changed.
webster_rr_parameters <- function(event) {
  target <- webster_event_targets(event)
  reference <- webster_population(ess_ratio = .5)
  u_source <- unname(reference$q[2])
  u_target <- unname(reference$e[2])
  ## Under this reference design both scores are Binomial(3,u).
  moment <- function(t, u) (1 - u + u * t)^3
  ratio <- function(t) moment(t, u_source) / moment(t, u_target)
  inverse_ratio <- function(r) {
    k <- r^(1/3)
    (1 - u_source - k * (1 - u_target)) / (k * u_target - u_source)
  }
  ## One free shape parameter remains after matching the six mean risks.
  ## .85 gives three prognostic EMs and valid cell probabilities in both events.
  baseline_rr_E <- .85
  baseline_rr_P <- inverse_ratio(target["AC", "p_C"] / target["BC", "p_C"] /
                                  ratio(baseline_rr_E))
  active_E <- vapply(c("p_A", "p_B"), function(arm)
    inverse_ratio(target["AC", arm] / target["BC", arm] / ratio(baseline_rr_P)), numeric(1))
  baseline_risk <- target["BC", "p_C"] /
    (moment(baseline_rr_E, u_target) * moment(baseline_rr_P, u_target))
  rr <- target["BC", c("p_A", "p_B")] /
    (baseline_risk * vapply(active_E, moment, numeric(1), u = u_target) *
       moment(baseline_rr_P, u_target))
  list(baseline_risk = unname(baseline_risk), rr_AC = unname(rr[1]), rr_BC = unname(rr[2]),
       em_rr_AC = unname(active_E[1] / baseline_rr_E),
       em_rr_BC = unname(active_E[2] / baseline_rr_E),
       baseline_rr_E = baseline_rr_E, baseline_rr_P = unname(baseline_rr_P))
}

webster_rr_design <- function(event = "common", baseline_risk = NULL,
                              rr_AC = NULL, rr_BC = NULL,
                              em_rr_AC = NULL, em_rr_BC = NULL,
                              baseline_rr_E = NULL, baseline_rr_P = NULL,
                              q = NULL, ess_ratio = .5,
                              direction = c(-1, 1, -1, 1, -1, 1)) {
  event <- webster_event(event)
  parameters <- list(baseline_risk = baseline_risk,
    rr_AC = rr_AC, rr_BC = rr_BC, em_rr_AC = em_rr_AC, em_rr_BC = em_rr_BC,
    baseline_rr_E = baseline_rr_E, baseline_rr_P = baseline_rr_P)
  supplied <- !vapply(parameters, is.null, logical(1))
  preset <- webster_rr_parameters(event)
  preset[supplied] <- parameters[supplied]
  parameters <- preset
  valid <- vapply(parameters, function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1))
  if (!all(valid) || parameters$baseline_risk > 1)
    stop("Risk multipliers must be finite positive scalars; baseline_risk must be in (0,1].")
  baseline_risk <- parameters$baseline_risk
  rr_AC <- parameters$rr_AC; rr_BC <- parameters$rr_BC
  em_rr_AC <- parameters$em_rr_AC; em_rr_BC <- parameters$em_rr_BC
  baseline_rr_E <- parameters$baseline_rr_E; baseline_rr_P <- parameters$baseline_rr_P
  population <- webster_population(q, ess_ratio, direction)
  e <- population$e; q <- population$q
  cols <- paste0("z", 1:6)
  tab <- expand.grid(S_E = 0:3, S_P = 0:3)
  tab$h_E <- (tab$S_E - .9) / sqrt(.63)
  tab$h_P <- (tab$S_P - .9) / sqrt(.63)
  tab$log_p_C <- log(baseline_risk) + log(baseline_rr_E) * tab$S_E +
    log(baseline_rr_P) * tab$S_P
  tab$log_RR_AC <- log(rr_AC) + log(em_rr_AC) * tab$S_E
  tab$log_RR_BC <- log(rr_BC) + log(em_rr_BC) * tab$S_E
  tab$p_C <- exp(tab$log_p_C)
  tab$p_A <- exp(tab$log_p_C + tab$log_RR_AC)
  tab$p_B <- exp(tab$log_p_C + tab$log_RR_BC)
  risks <- as.matrix(tab[c("p_C", "p_A", "p_B")])
  if (any(!is.finite(risks)) || any(risks <= 0 | risks > 1))
    stop("Invalid log-risk parameters: every C/A/B risk must be in (0,1]. No clipping is performed.")
  coefficient_truth <- rbind(
    AC = c(alpha = log(rr_AC) + .9 * log(em_rr_AC), lambda = sqrt(.63) * log(em_rr_AC)),
    BC = c(alpha = log(rr_BC) + .9 * log(em_rr_BC), lambda = sqrt(.63) * log(em_rr_BC)))
  list(dgm = "webster2024_rr_3em3nonem", effect_scale = "RR", parameters = parameters,
       e = e, q = q, ess_ratio = rel_ess_binary(q, e), event = event,
       score_center = .9, score_scale = sqrt(.63), risk_table = tab,
       marginal_truth = rbind(AC = rr_gop_marginal_truth(tab, q),
                             BC = rr_gop_marginal_truth(tab, e)),
       coefficient_truth = coefficient_truth,
       calibration = webster_calibration(event, tab, names(parameters)[supplied]),
       covariate_roles = list(effect_modifiers = cols[1:3],
                              non_effect_modifiers = cols[4:6]),
       source = list(doi = "10.1093/aje/kwae048",
         revision = "6d495d97b5002d16e5ae50bfd2f33c1a58f3dddc",
         file = "Sensitivity_analyses/simulating_variable_selection_RR.sas",
         adaptation = "Three EMs plus three prognostic non-EMs; two trials; original project covariates/ESS; balanced allocation."))
}
