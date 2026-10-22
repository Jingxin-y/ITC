#!/bin/bash
# The user's full 144-cell design; identical DGP/seed settings generate shared data.
set -euo pipefail
mode="${1:-dry-run}"
[[ "$#" -le 1 && ( "$mode" == submit || "$mode" == dry-run ) ]] || { echo "Usage: bash $0 [dry-run|submit]" >&2; exit 1; }
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
code_dir="${MAIC_CODE_DIR:-$script_dir}"
count=0
submit_cell() {
    local source="$1" level="$2" event="$3" weight="$4" n="$5" ess="$6"
    local dry=0
    [[ "$mode" == dry-run ]] && dry=1
    env MAIC_DRY_RUN="$dry" MAIC_DGP=ess_em MAIC_CODE_DIR="$code_dir" \
        MAIC_ESS_SOURCE="$source" MAIC_EM_LEVEL="$level" MAIC_EVENT="$event" MAIC_WEIGHT_MODE="$weight" \
        MAIC_EST_MODE=scores MAIC_EST_EM=z1,z2,z3 MAIC_EST_NON_EM=z4,z5,z6 MAIC_MODEL_WEIGHTING=unit \
        bash "$code_dir/submit_ess_em_grid.sh" "$n" "$ess"
    count=$((count+1))
}
for event in common rare; do
    for level in low high; do
        for n in 100 200 500; do
            for ess in 0.15 0.25 0.5 0.8; do submit_cell diffuse "$level" "$event" all "$n" "$ess"; done
        done
        for weight in em nonem; do
            for n in 100 500; do
                for ess in 0.15 0.5 0.8; do submit_cell diffuse "$level" "$event" "$weight" "$n" "$ess"; done
            done
        done
    done
    for source in em_driven nonem_driven; do
        for weight in all em nonem; do
            for n in 100 500; do
                for ess in 0.15 0.8; do submit_cell "$source" high "$event" "$weight" "$n" "$ess"; done
            done
        done
    done
done
echo "Total $mode scenarios: $count"
[[ "$count" == 144 ]] || { echo "Unexpected matrix size." >&2; exit 1; }
