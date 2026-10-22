#!/bin/bash
# One job = one ESS-source / EM-strength / calibration / event / n / ESS setting.
# RR risk targets v2 (2026-10-15); use this package's code snapshot.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${1:-}" == worker ]]; then
    : "${SLURM_JOB_ID:?Workers must be submitted through sbatch.}"
    code_dir="$2"; run_dir="$3"; n="$4"; ess="$5"; ci_mode="$6"; cores="$7"
    mode=worker
elif [[ "${1:-}" == collect ]]; then
    run_dir="${2:?Supply the result directory.}"
    code_dir="$run_dir/code_snapshot"; n=0; ess=0; ci_mode=full; cores=1
    mode=collect
else
    [[ "$#" == 2 ]] || { echo "Usage: bash $0 N ESS_RATIO" >&2; exit 1; }
    n="$1"; ess="$2"
    case "$n" in 100|200|500) ;; *) echo "N must be 100, 200 or 500." >&2; exit 1 ;; esac
    case "$ess" in
        0.15) ess_code=015 ;; 0.25) ess_code=025 ;;
        0.5) ess_code=050 ;; 0.8) ess_code=080 ;;
        *) echo "ESS_RATIO must be 0.15, 0.25, 0.5 or 0.8." >&2; exit 1 ;;
    esac
    code_dir="${MAIC_CODE_DIR:-$script_dir}"
    nsim="${MAIC_NSIM:-500}"; first_seed="${MAIC_FIRST_SEED:-1}"
    cores="${MAIC_CORES_PER_JOB:-96}"; mem_gb="${MAIC_MEM_GB:-20}"
    walltime="${MAIC_WALLTIME:-3:00:00}"; account="${MAIC_ACCOUNT:-def-liteep}"
    event="${MAIC_EVENT:-common}"; event="${event,,}"
    ess_source="${MAIC_ESS_SOURCE:-diffuse}"; ess_source="${ess_source,,}"
    em_level="${MAIC_EM_LEVEL:-high}"; em_level="${em_level,,}"
    weight_mode="${MAIC_WEIGHT_MODE:-all}"; weight_mode="${weight_mode,,}"
    model_weighting="${MAIC_MODEL_WEIGHTING:-unit}"; model_weighting="${model_weighting,,}"
    estimation_mode="${MAIC_EST_MODE:-scores}"; estimation_mode="${estimation_mode,,}"
    ci_mode="${MAIC_CI_MODE:-full}"
    [[ "${MAIC_DGP:-ess_em}" == ess_em ]] || { echo "This script requires MAIC_DGP=ess_em." >&2; exit 1; }
    for number in "$nsim" "$first_seed" "$cores" "$mem_gb"; do
        [[ "$number" =~ ^[1-9][0-9]*$ ]] || { echo "Counts, cores, memory and seeds must be positive integers." >&2; exit 1; }
    done
    [[ "$event" == common || "$event" == rare ]] || { echo "Invalid MAIC_EVENT." >&2; exit 1; }
    [[ "$ess_source" == diffuse || "$ess_source" == em_driven || "$ess_source" == nonem_driven ]] || { echo "Invalid MAIC_ESS_SOURCE." >&2; exit 1; }
    [[ "$em_level" == low || "$em_level" == high ]] || { echo "Invalid MAIC_EM_LEVEL." >&2; exit 1; }
    [[ "$weight_mode" == all || "$weight_mode" == em || "$weight_mode" == nonem ]] || { echo "Invalid MAIC_WEIGHT_MODE." >&2; exit 1; }
    [[ "$model_weighting" == unit || "$model_weighting" == maic ]] || { echo "Invalid MAIC_MODEL_WEIGHTING." >&2; exit 1; }
    [[ "$estimation_mode" == scores || "$estimation_mode" == none ]] || { echo "Invalid MAIC_EST_MODE." >&2; exit 1; }
    [[ "$ci_mode" == full || "$ci_mode" == smoke ]] || { echo "Invalid MAIC_CI_MODE." >&2; exit 1; }
    for file in load_anchored.R data/generate_webster_RR_scenarios.R estimation/cluster_all_methods.R estimation/cluster_grid.R estimation/helpers/CI_marginal_target.R; do
        [[ -f "$code_dir/$file" ]] || { echo "Missing code: $code_dir/$file" >&2; exit 1; }
    done
    export MAIC_DGP=ess_em MAIC_EVENT="$event" MAIC_ESS_SOURCE="$ess_source" MAIC_EM_LEVEL="$em_level"
    export MAIC_WEIGHT_MODE="$weight_mode" MAIC_MODEL_WEIGHTING="$model_weighting" MAIC_EST_MODE="$estimation_mode"
    export MAIC_EST_EM="${MAIC_EST_EM:-z1,z2,z3}" MAIC_EST_NON_EM="${MAIC_EST_NON_EM:-z4,z5,z6}"
    if [[ "${MAIC_DRY_RUN:-0}" == 1 ]]; then
        printf 'DRY RUN: event=%s source=%s EM=%s weights=%s n=%s full_ESS=%s nsim=%s CI=%s model=%s cores=%s time=%s\n' \
            "$event" "$ess_source" "$em_level" "$weight_mode" "$n" "$ess" "$nsim" "$ci_mode" "$model_weighting" "$cores" "$walltime"
        exit 0
    fi
    command -v sbatch >/dev/null || { echo "Run this on the cluster: sbatch is required." >&2; exit 1; }
    : "${SCRATCH:?Set SCRATCH to your cluster scratch directory.}"
    scenario="${event}_${ess_source}_${em_level}_weight_${weight_mode}_est_${estimation_mode}_fit_${model_weighting}_n${n}_ess${ess_code}"
    run_dir="$SCRATCH/MAIC_ess_em_v2/${scenario}_${ci_mode}_$(date +%Y%m%d_%H%M%S)_$$"
    mkdir -p "$run_dir/logs" "$run_dir/by_seed" "$run_dir/code_snapshot"
    cp -a -- "$code_dir/data" "$code_dir/estimation" "$code_dir/summary" "$code_dir/load_anchored.R" "$run_dir/code_snapshot/"
    cp -- "${BASH_SOURCE[0]}" "$run_dir/submit_ess_em_grid.sh"
    printf 'nsim=%s\nfirst_seed=%s\nci_mode=%s\nevent=%s\nestimation_mode=%s\nn=%s\ness_ratio=%s\ncores=%s\ndgp_mode=ess_em\ness_source=%s\nem_level=%s\nweight_mode=%s\nmodel_weighting=%s\n' \
        "$nsim" "$first_seed" "$ci_mode" "$event" "$estimation_mode" "$n" "$ess" "$cores" "$ess_source" "$em_level" "$weight_mode" "$model_weighting" > "$run_dir/submission.txt"
    job=$(sbatch --parsable --account="$account" --job-name="MAIC_${event}_${ess_source}_${em_level}_${weight_mode}_n${n}_ess${ess_code}" \
        --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="${mem_gb}G" --time="$walltime" \
        --export=ALL --output="$run_dir/logs/%j.out" --error="$run_dir/logs/%j.err" \
        "$run_dir/submit_ess_em_grid.sh" worker "$run_dir/code_snapshot" "$run_dir" "$n" "$ess" "$ci_mode" "$cores")
    job="${job%%;*}"; printf '%s\n' "$job" > "$run_dir/job_id.txt"
    echo "Submitted job $job: $scenario, $nsim replicates, $cores cores, time $walltime."
    echo "Results: $run_dir"
    exit 0
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
module load StdEnv/2023
module load gcc/12.3 "${MAIC_R_MODULE:-r/4.3.1}"
Rscript --vanilla - "$mode" "$code_dir" "$run_dir" "$n" "$ess" "$ci_mode" "$cores" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
mode <- args[1L]; code_dir <- args[2L]; run_dir <- args[3L]
source(file.path(code_dir,"estimation/cluster_all_methods.R"),encoding="UTF-8")
source(file.path(code_dir,"estimation/cluster_grid.R"),encoding="UTF-8")
if(mode=="collect") {
  maic_collect(run_dir)
} else {
  source(file.path(code_dir,"load_anchored.R"),encoding="UTF-8")
  source_anchored(code_dir)
  maic_run_scenario(code_dir,run_dir,as.integer(args[4L]),as.numeric(args[5L]),args[6L],as.integer(args[7L]))
  maic_collect(run_dir)
}
RSCRIPT
