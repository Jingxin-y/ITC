#!/bin/bash
set -euo pipefail
source "$RRNMW_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RRNMW_HOME"
sha256sum --quiet -c SHA256SUMS
exec 9>"$RRNMW_OUT/.locks/task_${SLURM_ARRAY_TASK_ID}.lock"
flock -n 9 || { echo 'This formal RR task already has an active worker'; exit 1; }
Rscript --vanilla run.R "$RRNMW_HOME" "$RRNMW_OUT" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
