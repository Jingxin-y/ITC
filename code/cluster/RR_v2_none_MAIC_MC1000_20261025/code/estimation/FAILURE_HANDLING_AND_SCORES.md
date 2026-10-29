# Failure handling and score standardization

The comparison table separates point-estimate `status` from `se_status`,
`interval_status`, and `p_status`. An `available` status with valid numerical
values contributes to the corresponding summary. Intervals also accept
`BC_boundary_truncated`, with a separate `n_boundary_truncated` count. Raw fit objects and BC diagnostic
endpoints remain available in the full comparison result for investigation.

- Unconverged BRM/FC/GLM fits do not contribute estimates or intervals to the main summary.
  Ordinary BRM also checks the final inner optimizer codes, rather than parameter
  movement alone. This is a convergence check, not a proof of a global optimum.
  Their final coefficients are retained in `nonconverged_estimates`, and raw fit
  objects in `nonconverged_fits`. BRM/FC objects include `last_optimizer` and
  `optimization_input` (data, weights, design matrices and final parameters as starts).
  Failed BRM/FC fits return these flagged objects without computing SEs or BC intervals.
- `fit_control` defaults to `brm_max_step=200`, `brm_optim_maxit=200`,
  `fc_max_step=300`, `glm_maxit=200`. One ordinary BRM outer step optimizes alpha
  with beta fixed, then beta with alpha fixed; each `optim()` call gets the inner
  limit of 200 and explicit `reltol=1e-6` (configurable as `brm_optim_reltol`).
  The actual inner tolerance is recorded in `optimization$optim_reltol` and
  nonconverged fits' `optimization_input$optim.reltol` for reproducible retries.
  FC instead solves the joint adjusted score using Newton steps.
  Default marginal BC uses `optim.maxit=200`; explicit `bc_target="coefficient"` retains `max.step=100`, `optim.maxit=300`.
  These are caps; each solver can stop earlier under its existing convergence rule.
  `brm_thres=1e-6` and `fc_thres=1e-6` are configurable separately. For a
  sensitivity run, use `fit_control=list(brm_thres=1e-3,fc_thres=1e-3)` in
  `run.anchored()`, `simulate_anchored()` or `maic_run_seed()`. BRM checks a
  squared relative parameter change AND both inner optimizer codes; FC checks
  its scaled adjusted score. These tolerances have different meanings. They
  do not change the inner optim tolerance, Exact score tolerance, root tolerance,
  or interval endpoint tolerance. The outer thresholds retain their defaults;
  the ordinary BRM inner objective tolerance is now explicitly `1e-6`.
- MAIC checks the optimizer code and final balance. Balance tolerance is `1e-6`
  per variable after dividing each error by `max(1, observed range)`.
  A failed calibration raises `maic_calibration_error` with diagnostics; the
  comparison catches it, keeps unweighted BRM, and marks weighted methods failed.
  Individual support checks are necessary but not sufficient; final residuals
  also catch infeasible joint targets. Low ESS alone is not a failure criterion.
- A BC interval truncated only because a confirmed rejected point was not found
  within the search range contributes its returned endpoints to interval counts
  and coverage, while keeping `BC_boundary_truncated` status. Failed searches,
  failed internal fits, invalid bounds, or missing diagnostics remain unavailable.
  A real failure at either endpoint takes precedence over truncation. P-value validity
  is tracked separately. The main pipeline now defaults to exact_marginal_target(); coefficient diagnostics below remain available in explicit coefficient mode.
- Failure rows retain expected parameter names. A method with no successful fits
  still appears in the summary with zero usable estimates and NA bias/coverage.
  Unknown truth keys remain errors. Coverage is conditional on available intervals;
  inspect its denominator and the parameter status table alongside coverage.
- The cluster marginal Wald exporter preserves a GLM marginal point separately
  when only variance computation fails. In that case `estimate/risk_C/risk_A`
  are retained, `status=se_status=p_status="se_failed"`, the original error is
  saved in `message`, and SE/CI stay NA. Point summaries may use a finite point
  with this status, but SE and coverage must apply their own validity filters.
  The marginal Wald exporter uses `MASS::ginv()` for both Ag and H, with its
  default tolerance `sqrt(.Machine$double.eps)`, without a QR rank veto.
  Singular matrices therefore do not automatically abort the calculation.
  Joint calculations record the inverse method, tolerance, reciprocal condition
  numbers and scaled Moore-Penrose residuals `A %*% Aplus %*% A - A` in `message`.
  These diagnostics do not change availability. The generalized inverse drops
  directions below its singular-value threshold; returning a finite SE does not
  establish identifiability or valid coverage. Failed decompositions, nonfinite
  inverses and nonpositive or nonfinite variances remain unavailable.
  Existing saved results are not rewritten.
  In particular, a singular joint bread does not imply an infinite variance.
  NA/NaN and numerical Inf must be counted separately from finite plot values
  outside the display range; neither is automatically a valid unbounded CI.
- Marginal Exact was restored on 2026-09-27 to the solver used for the
  2026-09-22 results: bounded L-BFGS-B, `alpha.bound=8`, `beta.bound=10`,
  `optim.reltol=1e-10`, `score.tol=1e-4`, and `optim.maxit=200`.
  Starts outside the nuisance bounds fail as in the original solver;
  `expand.bounds` is no longer accepted. Convergence checks the projected
  nuisance score and the marginal constraint. The original diagnostic fields
  are retained; see `helpers/CI_MARGINAL_EXACT.md`.
  The endpoint search still uses `psi.range`.
  Valid `boundary_truncated` intervals continue to count, while actual fitting
  failures remain failures. Explicit legacy coefficient/linear-target routines
  retain their original optimizers.
- `simulate_anchored(checkpoint_file=...)` saves after each completed replicate and
  preserves completed work if a later replicate raises an unexpected error. Such
  errors stop the run and are recorded with the interrupted seed. The checkpoint
  contains `seeds`, `requested_seeds`, and `complete`; automatic resume is not implemented.
  The CLI enables checkpoints by default and saves raw RDS/CSV before summarizing.
  It also writes `anchored_joint_nonconverged.csv` (last parameter values) and
  `anchored_joint_nonconverged.rds` (table plus raw fits and configuration), even
  when the table is empty. Main checkpoints retain the same candidates after each replicate.
  `anchored_joint_parameter_status.csv` records failures separately for each output.

## BC aggregate estimates versus score units

`get.aggre()` already uses sample means and sample covariance from the BC trial.
These are the AD moments used by MAIC and the joint sandwich, not population moments
from the generator. Custom-score AD moments are affine transformations of these
sample moments in `brm_score_summary()`.

Generation retains center 0.9 and scale sqrt(0.63) from the BC population truth.
Analysis defaults to `score_standardization="bc_estimated"` in the comparison,
single-replicate and batch entry points. `brm_scores_from_bc()` replaces the score
template's center/scale using reported BC sample moments, before any method is fit.
For a standalone fit, pass a resolved spec from `brm_scores_from_bc(spec,data.AD)`
or construct it with `make_brm_score_spec(...,data.ad=data.AD)`.
The explicit `score_standardization="fixed"` option retains prespecified units.

For an affine raw score S = b + a'Z, the analysis uses

    estimated_mean(S) = b + a' estimated_mean(Z)
    estimated_variance(S) = a' estimated_covariance(Z) a
    h = (S - estimated_mean(S)) / estimated_sd(S)

Use the same center and scale for the AC individuals and BC summaries. Marginal
means/variances alone do not identify the variance of a sum unless cross-covariances
are supplied or an explicit correlation assumption is made. Plugging sample
proportions into `reference_prob` still uses the constructor's independence formula;
it is not a substitute for the full sample covariance.

When center and scale are estimated afresh in each simulation replicate, coefficients
in those units also change. If the generating log RR is alpha + lambda*(S-c0)/s0,
then in h=(S-c_hat)/s_hat units it is

    alpha_star = alpha + lambda*(c_hat-c0)/s0
    lambda_star = lambda*s_hat/s0.

`reference_truth` stores these replicate-specific coefficient truths, and
`score_references` stores each realized center and scale. When a selected score cannot
represent the true linear RR, no exact coefficient truth is invented. With `truth=NULL`,
the summary uses available reference truths; custom truths can include a replicate `r`.
Bias, coverage and error SD compare each estimate to its corresponding truth.
`empirical_sd_estimate` separately reports the SD of estimates themselves, while
`empirical_sd` and SE ratios use error SD when truth is available. These summaries
target effects defined in each reported BC score scale, not a single fixed-scale coefficient.

Inference for a fixed population-standardized coefficient would need to propagate
the uncertainty in estimated score units, or transform inference back to fixed units.
The current joint sandwich accounts for calibration target means, but does not
automatically account for estimating a score's center/scale or their covariance with
the fit. BC still holds the score design and weights fixed. These are separate
extensions; this change does not claim that conditional BC intervals now propagate
all uncertainty in the two trials. Coverage validation remains necessary.
