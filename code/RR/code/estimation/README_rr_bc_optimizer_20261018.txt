RR BC numerical rescue update — 2026-10-18

Purpose
  Add the deployed RD deterministic multistart and nlminb rescue to marginal
  log-RR BC. The original initial L-BFGS-B attempt and two resume retries stay
  first and use their existing settings. Only after all original attempts fail
  are the new paths eligible. An originally successful fit is returned directly.

Rescue order (same as RD)
  1. Deterministic distinct starts: successful same-response neighboring target,
     nearby observed-data profile, cell-means start, and neutral start, when
     available. At most three additional starts are used by default.
  2. If none succeeds, bounded stats::nlminb uses the same likelihood and analytic
     gradient, starting from the best failed result and available alternatives.
  Cell means use half-event smoothing only to construct finite initial values.
  The responses, likelihood and derivatives do not receive added events.

Default optional controls accepted by the main bc_control adapter
  optim.multistart = TRUE
  optim.multistart.max = 3L
  optim.fallback = TRUE
  Set optim.multistart=FALSE and optim.fallback=FALSE to disable the new paths.
  Original retries and all original statistical/numerical thresholds remain.

Unchanged
  DGM, RR/OP probability model, ordinary likelihood, observed fixed anchor,
  candidate/anchor bootstrap statistic, MAIC fitting/standardization weights,
  bootstrap counts and RNG, parameter bounds, gradient/constraint acceptance,
  interval search and boundary handling, time limit, point/FC/Wald estimators.
  New retries are deterministic and reuse the same response and candidate effect.
  Permanently failed bootstrap fits remain failures; no draws are dropped or redrawn.

Diagnostics
  optimizer.version = rr_bc_multistart_nlminb_v1
  Attempt diagnostics include method, start label/vector, objective, projected
  score, constraint error and elapsed time. Recovery counters distinguish
  original resumes, multistart and nlminb. Source-file hashes change as expected.

Deployment
  Required runtime files:
    estimation/helpers/CI_marginal_target.R
    estimation/helpers/marginal_comparison.R
  Tests and this note are optional for cluster execution. Apply only to a code
  tree matching the original inputs backed up with this patch. Use a fresh run
  directory for updated simulations; do not mix old/new source hashes or assume
  existing checkpoints will recompute old BC results after replacing a source.
  Old result files and shared FC sources are not changed by this update.

Validation commands, run from the code root
  Rscript --vanilla estimation/tests/check_marginal_optimizer_rescue.R
  Rscript --vanilla estimation/tests/check_marginal_retry.R
  Rscript --vanilla estimation/tests/check_marginal_exact.R
  Rscript --vanilla estimation/tests/check_marginal_pipeline.R
  Rscript --vanilla estimation/tests/check_firth_solver_consistency.R
  The final check expects RD/ beside the RR source; an external RD module path
  can instead be supplied as its first argument. On Windows, configure LC_CTYPE
  as English_United States.utf8 if R cannot source the original UTF-8 comments.

Firth audit
  RR and standalone RD use the same max.likelihood.firth function body and
  firth_score_adjustment equation, mean-one weight normalization, total step
  budget (300 by default), score acceptance tolerance (1e-6 by default), Newton /
  SVD / line search and deterministic restarts. They deliberately use different
  RR versus RD probability mappings and first/second derivatives; the RR effect
  link is log(RR), whereas RD is atanh(RD). No Firth source was modified.

Scope of evidence
  Local regression, controlled failure/recovery and small bootstrap integration
  checks establish the numerical contract. They do not establish improved
  coverage or a particular return-rate gain in the full cluster simulation.
