#!/bin/bash
set -euo pipefail
source "$RDNM_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RDNM_HOME"
sha256sum --quiet -c SHA256SUMS
exec 9>"$RDNM_OUT/.locks/task_${SLURM_ARRAY_TASK_ID}.lock"
flock -n 9 || { echo 'This formal RD task already has an active worker'; exit 1; }
Rscript --vanilla run.R "$RDNM_HOME" "$RDNM_OUT" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
