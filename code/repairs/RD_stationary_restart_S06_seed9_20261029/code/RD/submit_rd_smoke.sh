#!/bin/bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RD_NSIM="${RD_NSIM:-2}" RD_CORES="${RD_CORES:-2}" RD_CI_MODE=smoke RD_WALLTIME="${RD_WALLTIME:-0:30:00}"
bash "$script_dir/submit_rd_grid.sh" "${1:-dry-run}" "${2:-100}" "${3:-0.5}"
