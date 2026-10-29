#!/bin/bash
set -euo pipefail
source "$RDFULL_HOME/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RDFULL_HOME"
sha256sum --quiet -c SHA256SUMS
exec 9>"$RDFULL_OUT/.locks/task_${SLURM_ARRAY_TASK_ID}.lock"
flock -n 9 || { echo 'This formal RD task already has an active worker'; exit 1; }
Rscript --vanilla run.R "$RDFULL_HOME" "$RDFULL_OUT" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
