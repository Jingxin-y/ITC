## Shared population definition and calibration targets for direct-risk DGMs.
webster_event <- function(event) {
  if (!is.character(event) || length(event) != 1L || is.na(event) ||
      !tolower(event) %in% c("common", "rare"))
    stop('event must be "common" or "rare" (case-insensitive).')
  tolower(event)
}

webster_population <- function(q = NULL, ess_ratio = .5,
                               direction = c(-1, 1, -1, 1, -1, 1)) {
  e <- c(.7, .3, .7, .3, .7, .3)
  if (is.null(q)) {
    if (!is.numeric(ess_ratio) || length(ess_ratio) != 1L ||
        !is.finite(ess_ratio) || ess_ratio <= 0 || ess_ratio > 1)
      stop("ess_ratio must be a finite scalar in (0,1].")
    if (!is.numeric(direction) || length(direction) != 6L ||
        any(!is.finite(direction)) || !any(direction != 0))
      stop("direction must contain six finite numbers and a nonzero direction.")
    q <- ipd_prob_for_ess_theory(ess_ratio, e, direction)
  }
  cols <- paste0("z", 1:6)
  if (!is.numeric(q) || length(q) != 6L || any(!is.finite(q)) || any(q <= 0 | q >= 1))
    stop("q must contain six finite probabilities strictly between 0 and 1.")
  if (!is.null(names(q))) {
    if (anyDuplicated(names(q)) || !setequal(names(q), cols))
      stop("q must be unnamed or named z1,...,z6.")
    q <- q[cols]
  }
  names(e) <- names(q) <- cols
  list(e = e, q = q, ess_ratio = rel_ess_binary(q, e))
}

## Exact old RR-GOP marginal risks, recorded at ESS=.5/default direction.
## alpha_AC=.20, alpha_BC=.19, lambda=.30, rho=.50,
## beta_E=.20, beta_P=.40, beta0=-.26 (common) / -8.53 (rare).
## These constants are calibration targets, not an OP/GOP outcome generator.
webster_event_targets <- function(event) {
  event <- webster_event(event)
  values <- if (event == "common") c(
    .400992802107635, .418596226792393,
    .599045919639806, .514669087904676,
    .529272607280753, .502197586315061) else c(
    .0501312948679940, .0488309131792698,
    .0762860380427614, .0608090036668075,
    .0667672466456628, .0589376565366341)
  matrix(values, nrow = 2L,
         dimnames = list(c("AC", "BC"), c("p_C", "p_A", "p_B")))
}

webster_calibration <- function(event, tab, overridden_parameters) {
  reference <- webster_population(ess_ratio = .5)
  reference_risks <- rbind(AC = rr_gop_marginal_truth(tab, reference$q),
                           BC = rr_gop_marginal_truth(tab, reference$e))[, c("p_C", "p_A", "p_B")]
  targets <- webster_event_targets(event)
  list(reference_ess_ratio = .5, reference_direction = c(-1, 1, -1, 1, -1, 1),
       reference_q = reference$q, reference_e = reference$e,
       target_risks = targets, reference_risks = reference_risks,
       matches_reference_targets = max(abs(reference_risks - targets)) < 1e-11,
       overridden_parameters = overridden_parameters,
       scope = "Fixed presets calibrated to original AC and BC populations at ESS=.5; changing ESS/q does not recalibrate the risk table.")
}

webster_rd_marginal_truth <- function(tab, prob) {
  truth <- rr_gop_marginal_truth(tab, prob)
  c(truth, RD_AC = unname(truth["p_A"] - truth["p_C"]),
    RD_BC = unname(truth["p_B"] - truth["p_C"]),
    RD_AB = unname(truth["p_A"] - truth["p_B"]))
}
