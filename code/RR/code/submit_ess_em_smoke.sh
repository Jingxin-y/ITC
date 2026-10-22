#!/bin/bash
# Six real small pipeline cases, including marginal BC intervals; one Slurm job.
set -euo pipefail
if [[ "${1:-}" == worker ]]; then
    : "${SLURM_JOB_ID:?Workers must be submitted through sbatch.}"
    code_dir="$2"; run_dir="$3"; cores="$4"
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
    module load StdEnv/2023
    module load gcc/12.3 "${MAIC_R_MODULE:-r/4.3.1}"
    Rscript --vanilla "$code_dir/estimation/tests/run_cluster_ess_em_smoke.R" "$code_dir" "$run_dir" "$cores"
    exit 0
fi
[[ "$#" == 0 ]] || { echo "Usage: bash $0" >&2; exit 1; }
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
code_dir="${MAIC_CODE_DIR:-$script_dir}"
cores="${MAIC_SMOKE_CORES:-2}"; mem_gb="${MAIC_SMOKE_MEM_GB:-8}"
walltime="${MAIC_SMOKE_WALLTIME:-1:00:00}"; account="${MAIC_ACCOUNT:-def-liteep}"
for number in "$cores" "$mem_gb"; do
    [[ "$number" =~ ^[1-9][0-9]*$ ]] || { echo "Smoke cores and memory must be positive integers." >&2; exit 1; }
done
command -v sbatch >/dev/null || { echo "Run this on the cluster: sbatch is required." >&2; exit 1; }
: "${SCRATCH:?Set SCRATCH to your cluster scratch directory.}"
[[ -f "$code_dir/estimation/tests/run_cluster_ess_em_smoke.R" ]] || { echo "Missing smoke test." >&2; exit 1; }
run_dir="$SCRATCH/MAIC_ess_em_v2/smoke_$(date +%Y%m%d_%H%M%S)_$$"
mkdir -p "$run_dir/logs" "$run_dir/code_snapshot"
cp -a -- "$code_dir/data" "$code_dir/estimation" "$code_dir/summary" "$code_dir/load_anchored.R" "$run_dir/code_snapshot/"
cp -- "${BASH_SOURCE[0]}" "$run_dir/submit_ess_em_smoke.sh"
job=$(sbatch --parsable --account="$account" --job-name=MAIC_ess_em_smoke \
    --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="${mem_gb}G" --time="$walltime" \
    --export=ALL --output="$run_dir/logs/%j.out" --error="$run_dir/logs/%j.err" \
    "$run_dir/submit_ess_em_smoke.sh" worker "$run_dir/code_snapshot" "$run_dir" "$cores")
job="${job%%;*}"; printf '%s\n' "$job" > "$run_dir/job_id.txt"
echo "Smoke job: $job"
echo "Results: $run_dir"
printf 'After completion:\n  cat %q\n  cat %q\n' "$run_dir/SMOKE_PASS.txt" "$run_dir/logs/$job.out"
