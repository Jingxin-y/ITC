## Independent ADDITIVE risk DGM inspired by the paper's RD setting.
## Project adaptation: three EMs, three prognostic non-EMs, AC/BC trials.
## It is a data generator only; the current RR comparison is not an RD analysis.
## p_C=b+beta_E*S_E+beta_P*S_P
## p_A=p_C+rd_AC+em_rd_AC*S_E; p_B=p_C+rd_BC+em_rd_BC*S_E
## Consequently S_P is not an EM on either the conditional or marginal RD scale.
webster_rd_parameters <- function(event) {
  event <- webster_event(event)
  target <- webster_event_targets(event)
  reference <- webster_population(ess_ratio = .5)
  mean_source <- sum(c(1-reference$q[1], reference$q[2], 1-reference$q[3]))
  mean_target <- .9
  delta_mean <- mean_source - mean_target
  ## One free baseline shape parameter; keep all three non-EMs prognostic.
  beta_P <- if (event == "common") .02 else .002
  beta_E <- (target["AC", "p_C"] - target["BC", "p_C"]) / delta_mean - beta_P
  b <- target["BC", "p_C"] - (beta_E + beta_P) * mean_target
  rd_source <- target["AC", c("p_A", "p_B")] - target["AC", "p_C"]
  rd_target <- target["BC", c("p_A", "p_B")] - target["BC", "p_C"]
  em <- (rd_source - rd_target) / delta_mean
  rd <- rd_target - em * mean_target
  list(baseline_risk = unname(b), rd_AC = unname(rd[1]), rd_BC = unname(rd[2]),
       em_rd_AC = unname(em[1]), em_rd_BC = unname(em[2]),
       beta_E = unname(beta_E), beta_P = beta_P)
}

webster_rd_design <- function(event = "common", baseline_risk = NULL,
                              rd_AC = NULL, rd_BC = NULL,
                              em_rd_AC = NULL, em_rd_BC = NULL,
                              beta_E = NULL, beta_P = NULL,
                              q = NULL, ess_ratio = .5,
                              direction = c(-1, 1, -1, 1, -1, 1)) {
  event <- webster_event(event)
  parameters <- list(baseline_risk = baseline_risk, rd_AC = rd_AC, rd_BC = rd_BC,
    em_rd_AC = em_rd_AC, em_rd_BC = em_rd_BC, beta_E = beta_E, beta_P = beta_P)
  supplied <- !vapply(parameters, is.null, logical(1))
  preset <- webster_rd_parameters(event)
  preset[supplied] <- parameters[supplied]
  parameters <- preset
  valid <- vapply(parameters, function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))
  if (!all(valid)) stop("Additive-risk coefficients must be finite numeric scalars.")
  population <- webster_population(q, ess_ratio, direction)
  e <- population$e; q <- population$q
  tab <- expand.grid(S_E = 0:3, S_P = 0:3)
  tab$h_E <- (tab$S_E - .9) / sqrt(.63)
  tab$h_P <- (tab$S_P - .9) / sqrt(.63)
  tab$p_C <- parameters$baseline_risk + parameters$beta_E * tab$S_E +
    parameters$beta_P * tab$S_P
  tab$RD_AC <- parameters$rd_AC + parameters$em_rd_AC * tab$S_E
  tab$RD_BC <- parameters$rd_BC + parameters$em_rd_BC * tab$S_E
  tab$p_A <- tab$p_C + tab$RD_AC
  tab$p_B <- tab$p_C + tab$RD_BC
  risks <- as.matrix(tab[c("p_C", "p_A", "p_B")])
  if (any(!is.finite(risks)) || any(risks <= 0 | risks >= 1))
    stop("Invalid additive-risk parameters: every C/A/B risk must be in (0,1). No clipping is performed.")
  ## Actual RR columns are descriptive; do not assume a linear log-RR truth.
  tab$log_RR_AC <- log(tab$p_A / tab$p_C)
  tab$log_RR_BC <- log(tab$p_B / tab$p_C)
  coefficient_truth <- rbind(
    AC = c(intercept = parameters$rd_AC + .9 * parameters$em_rd_AC,
           slope = sqrt(.63) * parameters$em_rd_AC),
    BC = c(intercept = parameters$rd_BC + .9 * parameters$em_rd_BC,
           slope = sqrt(.63) * parameters$em_rd_BC))
  cols <- paste0("z", 1:6)
  list(dgm = "webster2024_rd_3em3nonem", effect_scale = "RD", parameters = parameters,
       e = e, q = q, ess_ratio = population$ess_ratio, event = event,
       score_center = .9, score_scale = sqrt(.63), risk_table = tab,
       marginal_truth = rbind(AC = webster_rd_marginal_truth(tab, q),
                             BC = webster_rd_marginal_truth(tab, e)),
       rd_coefficient_truth = coefficient_truth,
       calibration = webster_calibration(event, tab, names(parameters)[supplied]),
       covariate_roles = list(effect_modifiers = cols[1:3], non_effect_modifiers = cols[4:6]),
       source = list(doi = "10.1093/aje/kwae048",
         adaptation = "Direct additive risk; three EMs plus three prognostic non-EMs; two trials; original populations/ESS and marginal-event targets. RD generation only; no change to RR estimation."))
}
