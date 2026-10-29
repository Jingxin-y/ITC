#!/bin/bash
set -euo pipefail
source "$RRFULL_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RRFULL_HOME"
sha256sum --quiet -c SHA256SUMS
exec 9>"$RRFULL_OUT/.locks/task_${SLURM_ARRAY_TASK_ID}.lock"
flock -n 9 || { echo 'This formal RR task already has an active worker'; exit 1; }
Rscript --vanilla run.R "$RRFULL_HOME" "$RRFULL_OUT" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
