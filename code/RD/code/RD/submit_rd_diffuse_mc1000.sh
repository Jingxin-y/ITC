#!/bin/bash
# Same MC1000 diffuse groups as submit_diffuse_mc1000.sh; RD methods, 15 GB.
set -euo pipefail

usage() {
    echo "Usage: bash $0 {all|em|nonem} {scores|none} [dry-run|submit]" >&2
}
[[ "$#" -ge 2 && "$#" -le 3 ]] || { usage; exit 1; }
weight_mode="$1"; estimation_mode="$2"; action="${3:-dry-run}"
[[ "$action" == dry-run || "$action" == submit ]] || { usage; exit 1; }

case "$weight_mode" in
    all|nonem) sizes=(100 500); ess_values=(0.15 0.5 0.8); expected=24 ;;
    em) sizes=(100 200 500); ess_values=(0.15 0.25 0.5 0.8); expected=48 ;;
    *) usage; exit 1 ;;
esac
case "$estimation_mode" in
    scores) model_weighting=unit; cores=32; walltime=12:00:00 ;;
    none) model_weighting=maic; cores=16; walltime=03:00:00 ;;
    *) usage; exit 1 ;;
esac

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
code_dir="${RD_CODE_DIR:-$(dirname -- "$script_dir")}";
runner="$script_dir/submit_rd_grid.sh"
[[ -f "$runner" ]] || { echo "Missing runner: $runner" >&2; exit 1; }
if [[ "$action" == submit ]]; then
    : "${SCRATCH:?Set SCRATCH to your cluster scratch directory.}"
fi
result_root="${RD_RESULT_ROOT:-${SCRATCH:-/SCRATCH}/MAIC_RD_MC1000_20261007}"
echo "RD group: calibration=$weight_mode model=$estimation_mode/$model_weighting action=$action"
echo "Results root: $result_root"

count=0
for event in rare common; do
    for level in low high; do
        for n in "${sizes[@]}"; do
            for ess in "${ess_values[@]}"; do
                if ! env RD_CODE_DIR="$code_dir" RD_RESULT_ROOT="$result_root" \
                    RD_ESS_SOURCE=diffuse RD_EM_LEVEL="$level" RD_EVENT="$event" \
                    RD_WEIGHT_MODE="$weight_mode" RD_EST_MODE="$estimation_mode" \
                    RD_MODEL_WEIGHTING="$model_weighting" \
                    RD_NSIM=1000 RD_FIRST_SEED=1 RD_CI_MODE=full \
                    RD_CORES="$cores" RD_MEM_GB=15 RD_WALLTIME="$walltime" \
                    bash "$runner" "$action" "$n" "$ess"; then
                    echo "Failed: $event EM=$level weights=$weight_mode model=$estimation_mode n=$n ESS=$ess; stopping." >&2
                    exit 1
                fi
                count=$((count+1))
            done
        done
    done
done
echo "Total RD $action scenarios: $count (expected $expected)"
[[ "$count" == "$expected" ]] || { echo "Unexpected scenario count." >&2; exit 1; }
