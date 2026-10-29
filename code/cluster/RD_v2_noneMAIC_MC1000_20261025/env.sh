#!/bin/bash
export RDNM_HOME
RDNM_HOME=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export RDNM_OUT="${RDNM_OUT:-${SCRATCH:-/scratch/$USER}/RD_v2_noneMAIC_MC1000_results_20261025}"
export RDNM_ACCOUNT="${RDNM_ACCOUNT:-def-liteep}"
# Optional overrides: RDNM_CPUS, RDNM_MEM, RDNM_TIME. Defaults: 64 / 10G / 03:00:00.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
