# RD BC optimizer update, 2026-10-12

This update changes only numerical constrained fitting inside RD BC. The DGM,
saturated atanh-RD/OP model, ordinary weighted likelihood, fixed observed anchor,
MAIC weights, probability mapping, point/FC/Wald estimators, bootstrap budgets,
parameter bounds and score/constraint acceptance tolerances remain unchanged.
The proposed natural-linear RD model and the risk-space prototype are not used.

## Solver sequence

1. Run the existing L-BFGS-B fit and its two resume retries.
2. On exhaustion, try up to three distinct deterministic starts: a successful
   neighboring target for the same response, a nearby observed-data profile,
   arm/cell means, and a neutral start, in that order when available.
   Compare the verified solutions and retain the lowest negative log likelihood.
3. If all those attempts fail, use bounded `stats::nlminb` with the same analytic
   gradient, starting from the best failed solution and available alternatives.
   Stop at the first verified fallback solution. No extra package is needed.

The arm/cell means use half-event smoothing ONLY to produce finite initial
parameters. The responses, counts, weights, objective and derivatives supplied
to every optimization attempt are unmodified. No event is added to a fit.
Missing arm/cell totals do not receive invented counts or a cell-means start.

An optimizer return code of zero alone is insufficient. Both algorithms must
satisfy the existing projected-score and marginal-constraint tests. Timeouts
propagate. Exhausted failures stay failed: bootstrap draws are never dropped,
replaced, or resampled. Neighbor storage is bounded to 16 observed-data starts
and one previous response; it does not retain an entire bootstrap history.

`nlminb` reference: [official R documentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlminb.html).

## Controls and diagnostics

Defaults are enabled through the existing `bc_control` interface:

```r
bc_control = list(optim.multistart = TRUE,
                  optim.multistart.max = 3L,
                  optim.fallback = TRUE)
```

Set both booleans to `FALSE` to compare with the previous numerical path.
`optim.multistart.max` permits integers 1 through 8. Retry iteration limits and
the time limit retain their existing controls. Full MC1000 need not be run to
test the update; start with a separate small batch of the difficult rare cases.

Diagnostics add `optimizer.version="rd_bc_multistart_nlminb_v1"`, attempt method,
start label, start vector, objective, projected score, constraint error and time.
Counts include `multistart.profile.fits`, `fallback.profile.fits`,
`multistart.attempts`, `fallback.attempts`, `recovered.by.multistart` and
`recovered.by.fallback`. Existing failed/retried response diagnostics remain.

## Cluster update

The only required replacement is:

```text
RD/estimation/CI_marginal_rd_target.R
```

The two test scripts and this note are optional. Keep the old solver outside
the `RD` source tree if you want to compare it, then upload the replacement.
Use a new result root. Checkpoint hashes intentionally reject results produced
by the prior implementation, even though the data-generating design is unchanged.
Existing submission scripts and their resource settings work without changes.

After uploading the solver and test script:

```bash
cd "$HOME/MAIC_ess_em_20261004"
Rscript --vanilla RD/tests/check_rd_optimizer_rescue.R "$PWD"
```

An optional paired numerical benchmark (five cases, n=100, ESS=.15, diffuse;
20 draws per candidate) can run as:

```bash
Rscript --vanilla RD/tests/compare_rd_optimizer_rescue.R \
  "$PWD" "$SCRATCH/MAIC_RD_BC_optimizer_probe_20261012" 20
```

If a saved pre-update solver exists at
`RD/validation/CI_marginal_rd_target_before.R`, the benchmark reads it; otherwise
it disables the new rescue paths as its baseline. It compares identical
candidate-specific bootstrap responses and reports every fit failure. It
separately records unavailable original point fits and observed candidate fits.
This is a numerical profile benchmark, not a coverage or full-CI return study.

## Validation limits

The forced-failure tests prove that distinct starts and nlminb can recover
numerical failures while preserving the objective, responses, RNG, bounds and
acceptance criteria. They also prove that permanent failures remain visible.
They do not prove an improved return rate for the original cluster batch.
The paired natural-data benchmark and integration results are documented in
the accompanying validation report. Zero/all-event boundary behavior and
nonidentifiable designs may still fail. Coverage needs a subsequent simulation.
