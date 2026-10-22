#!/bin/bash
# Same 144 cells as the existing RR ESS/EM matrix. Default is dry-run.
set -euo pipefail
mode="${1:-dry-run}"
[[ "$#" -le 1 && ( "$mode" == dry-run || "$mode" == submit ) ]] || exit 1
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
count=0
cell() {
    env RD_ESS_SOURCE="$1" RD_EM_LEVEL="$2" RD_EVENT="$3" RD_WEIGHT_MODE="$4" RD_MODEL_WEIGHTING=unit \
      bash "$script_dir/submit_rd_grid.sh" "$mode" "$5" "$6"
    count=$((count+1))
}
for event in common rare; do
  for level in low high; do
    for n in 100 200 500; do
      for ess in 0.15 0.25 0.5 0.8; do cell diffuse "$level" "$event" all "$n" "$ess"; done
    done
    for weight in em nonem; do
      for n in 100 500; do
        for ess in 0.15 0.5 0.8; do cell diffuse "$level" "$event" "$weight" "$n" "$ess"; done
      done
    done
  done
  for source in em_driven nonem_driven; do
    for weight in all em nonem; do
      for n in 100 500; do
        for ess in 0.15 0.8; do cell "$source" high "$event" "$weight" "$n" "$ess"; done
      done
    done
  done
done
echo "Total RD $mode scenarios: $count"
[[ "$count" == 144 ]]
