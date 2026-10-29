#!/bin/bash
# Default is read-only dry-run. Explicit submit creates a cluster job.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mode="${1:-dry-run}"
if [[ "$mode" == worker ]]; then
    : "${SLURM_JOB_ID:?Workers must run under sbatch.}"
    code_dir="$2"; rd_dir="$code_dir/RD"; run_dir="$3"
    n="$4"; ess="$5"; cores="$6"; nsim="$7"; first_seed="$8"
elif [[ "$mode" == collect ]]; then
    run_dir="${2:?Supply result directory}"; code_dir="$run_dir/code_snapshot"; rd_dir="$code_dir/RD"
    n=0; ess=0; cores=1; nsim=1; first_seed=1
else
    [[ "$mode" == dry-run || "$mode" == submit ]] || { echo "Usage: bash $0 [dry-run|submit|collect] [N] [ESS]" >&2; exit 1; }
    n="${2:-100}"; ess="${3:-0.5}"
    case "$n" in 100|200|500) ;; *) echo "N must be 100, 200 or 500" >&2; exit 1;; esac
    case "$ess" in 0.15|0.25|0.5|0.8|1) ;; *) echo "Unsupported ESS" >&2; exit 1;; esac
    code_dir="${RD_CODE_DIR:-$(dirname -- "$script_dir")}"; rd_dir="$script_dir"
    nsim="${RD_NSIM:-500}"; first_seed="${RD_FIRST_SEED:-1}"; cores="${RD_CORES:-32}"
    for v in "$nsim" "$first_seed" "$cores"; do [[ "$v" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid count/seed/cores" >&2; exit 1; }; done
    (( first_seed + nsim - 1 + 1000000 <= 2147483647 )) || { echo "Seed exceeds R integer range" >&2; exit 1; }
    export RD_EVENT="${RD_EVENT:-common}" RD_ESS_SOURCE="${RD_ESS_SOURCE:-diffuse}" RD_EM_LEVEL="${RD_EM_LEVEL:-high}"
    export RD_WEIGHT_MODE="${RD_WEIGHT_MODE:-all}" RD_MODEL_WEIGHTING="${RD_MODEL_WEIGHTING:-unit}" RD_CI_MODE="${RD_CI_MODE:-full}"
    export RD_EST_MODE="${RD_EST_MODE:-scores}"
    export RD_BRM_EFFECT_MODEL="${RD_BRM_EFFECT_MODEL:-linear_atanh}"
    [[ "$RD_EVENT" == common || "$RD_EVENT" == rare ]] || exit 1
    [[ "$RD_ESS_SOURCE" == diffuse || "$RD_ESS_SOURCE" == em_driven || "$RD_ESS_SOURCE" == nonem_driven ]] || exit 1
    [[ "$RD_EM_LEVEL" == low || "$RD_EM_LEVEL" == high ]] || exit 1
    [[ "$RD_WEIGHT_MODE" == all || "$RD_WEIGHT_MODE" == em || "$RD_WEIGHT_MODE" == nonem ]] || exit 1
    [[ "$RD_MODEL_WEIGHTING" == unit || "$RD_MODEL_WEIGHTING" == maic ]] || exit 1
    [[ "$RD_EST_MODE" == scores || "$RD_EST_MODE" == none ]] || { echo "RD_EST_MODE must be scores or none" >&2; exit 1; }
    [[ "$RD_BRM_EFFECT_MODEL" == linear_atanh || "$RD_BRM_EFFECT_MODEL" == saturated ]] || { echo "RD_BRM_EFFECT_MODEL must be linear_atanh or saturated" >&2; exit 1; }
    brm_tag="$RD_BRM_EFFECT_MODEL"
    [[ "$RD_EST_MODE" != none ]] || brm_tag=intercept_only
    [[ "$RD_CI_MODE" == smoke || "$RD_CI_MODE" == full || "$RD_CI_MODE" == wald ]] || exit 1
    for f in load_anchored.R data/generate_webster_RD.R data/generate_webster_RR_scenarios.R; do
        [[ -f "$code_dir/$f" ]] || { echo "Missing shared code: $f" >&2; exit 1; }
    done
    if [[ "$mode" == dry-run ]]; then
        printf 'RD DRY RUN event=%s source=%s strength=%s weights=%s est=%s model=%s brm=%s n=%s ESS=%s nsim=%s CI=%s cores=%s mem=%sG\n' \
          "$RD_EVENT" "$RD_ESS_SOURCE" "$RD_EM_LEVEL" "$RD_WEIGHT_MODE" "$RD_EST_MODE" "$RD_MODEL_WEIGHTING" "$brm_tag" "$n" "$ess" "$nsim" "$RD_CI_MODE" "$cores" "${RD_MEM_GB:-10}"
        exit 0
    fi
    command -v sbatch >/dev/null || { echo "sbatch required for submit" >&2; exit 1; }
    : "${SCRATCH:?Set SCRATCH to the cluster scratch directory}"
    scenario="RD_${RD_EVENT}_${RD_ESS_SOURCE}_${RD_EM_LEVEL}_weight_${RD_WEIGHT_MODE}_est_${RD_EST_MODE}_fit_${RD_MODEL_WEIGHTING}_brm_${brm_tag}_n${n}_ess${ess}"
    result_root="${RD_RESULT_ROOT:-$SCRATCH/MAIC_RD}"
    run_dir="$result_root/${scenario}_${RD_CI_MODE}_$(date +%Y%m%d_%H%M%S)_$$"
    mkdir -p "$run_dir/logs" "$run_dir/by_seed" "$run_dir/code_snapshot/RD"
    cp -a -- "$code_dir/data" "$code_dir/estimation" "$code_dir/summary" "$code_dir/load_anchored.R" "$run_dir/code_snapshot/"
    if [[ -f "$code_dir/BRMPLUS_LICENSE.txt" ]]; then cp -- "$code_dir/BRMPLUS_LICENSE.txt" "$run_dir/code_snapshot/"; fi
    cp -a -- "$rd_dir/data" "$rd_dir/estimation" "$rd_dir/summary" "$rd_dir/tests" "$run_dir/code_snapshot/RD/"
    cp -- "$rd_dir/"*.R "$rd_dir/"*.sh "$rd_dir/README.md" "$run_dir/code_snapshot/RD/"
    cp -- "${BASH_SOURCE[0]}" "$run_dir/submit_rd_grid.sh"
    printf 'effect_scale=RD\nscenario=%s\nnsim=%s\nfirst_seed=%s\nestimation_mode=%s\nbrm_effect_model=%s\ncores=%s\nmem_gb=%s\n' "$scenario" "$nsim" "$first_seed" "$RD_EST_MODE" "$brm_tag" "$cores" "${RD_MEM_GB:-10}" > "$run_dir/submission.txt"
    job=$(sbatch --parsable --account="${RD_ACCOUNT:-def-liteep}" --job-name="$scenario" \
      --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="${RD_MEM_GB:-10}G" --time="${RD_WALLTIME:-3:00:00}" \
      --export=ALL --output="$run_dir/logs/%j.out" --error="$run_dir/logs/%j.err" \
      "$run_dir/submit_rd_grid.sh" worker "$run_dir/code_snapshot" "$run_dir" "$n" "$ess" "$cores" "$nsim" "$first_seed")
    printf '%s\n' "$job" > "$run_dir/job_id.txt"
    echo "Submitted RD job $job; results: $run_dir"
    exit 0
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
module load StdEnv/2023
module load gcc/12.3 "${RD_R_MODULE:-r/4.3.1}"
Rscript --vanilla - "$mode" "$code_dir" "$rd_dir" "$run_dir" "$n" "$ess" "$cores" "$nsim" "$first_seed" <<'RSCRIPT'
a <- commandArgs(trailingOnly=TRUE)
source(file.path(a[3],"load_anchored_rd.R"),encoding="UTF-8")
source_anchored_rd(a[2],a[3])
if(a[1]=="collect") rd_collect(a[4]) else {
  seeds <- seq.int(as.integer(a[9]),length.out=as.integer(a[8]))
  rd_run_scenario(a[2],a[3],a[4],n=as.integer(a[5]),ess_ratio=as.numeric(a[6]),
    seeds=seeds,cores=as.integer(a[7]),event=Sys.getenv("RD_EVENT"),
    ess_source=Sys.getenv("RD_ESS_SOURCE"),em_level=Sys.getenv("RD_EM_LEVEL"),
    weight_mode=Sys.getenv("RD_WEIGHT_MODE"),model_weighting=Sys.getenv("RD_MODEL_WEIGHTING"),
    ci_mode=Sys.getenv("RD_CI_MODE"),estimation_mode=Sys.getenv("RD_EST_MODE","scores"),
    brm_effect_model=Sys.getenv("RD_BRM_EFFECT_MODEL","linear_atanh"))
}
RSCRIPT
