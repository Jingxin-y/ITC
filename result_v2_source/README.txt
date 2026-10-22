ITC RR / RD source sharing bundle — result_v2 version, revised 2026-10-22

RR/code and RD/code contain the complete, unchanged runtime sources.
CI_methods contains a separate review guide and exact source copies for
non-Wald intervals only: RR BC, RD BC and RD MN-ESS.
plotting contains the original figure helpers and portable saved-result entry points.
config/SCENARIOS.csv contains configuration parameters only.
metadata contains source origins, checksums and validation scope.

Historical numerical results, figures, audit reports, older snapshots and repair
archives have been removed from this sharing bundle. Wald code remains in the
complete runtime but is not separately listed in the CI review guide.

RR is an assembled copy of v2 DGM sources plus previously published BC/FC fixes.
RD preserves the 128 frozen runtime files. No new estimator or CI algorithm was
introduced, and full_CI_300_500 experiments are not included.

See README_zh.txt, CI_methods/README_zh.md and plotting/README_zh.txt.
Run python verify_checksums.py to verify the release.
The complete per-replicate CSV/RDS inputs needed for plotting are not included.
