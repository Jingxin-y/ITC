#!/bin/sh
#SBATCH --job-name=RLP_all_alt
#SBATCH --cpus-per-task=96
#SBATCH --time=2:00:00
#SBATCH --mem=20G

set -eu

export code_dir="${code_dir:-/home/yanjin41/brmplus_simulation/compareITC/sandwich}"
export result_dir="${result_dir:-/scratch/yanjin41/RRRDOR/brmplus_simulation/ITC/sandwich}"

R="${1:-1000}"
RUN_EXACT="${2:-TRUE}"
SKIP_EXISTING="${3:-TRUE}"
INCLUDE_GCOMP="${INCLUDE_GCOMP:-TRUE}"

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

mkdir -p "$result_dir/Rout"

run_file="$code_dir/run_simulation_anchored_RLP.R"
job_id="${SLURM_JOB_ID:-manual}"
rout_file="$result_dir/Rout/RR_simucpp_anchored_RLP_all_null_R_${R}_${job_id}.Rout"

if [ ! -f "$run_file" ]; then
  echo "Cannot find run file: $run_file" >&2
  exit 1
fi

r_args="R=${R} run_exact=${RUN_EXACT} skip_existing=${SKIP_EXISTING} include_gcomp=${INCLUDE_GCOMP} result_dir='${result_dir}' code_dir='${code_dir}'"

echo "Running all anchored RLP null scenarios"
echo "run_file=$run_file"
echo "result_dir=$result_dir"
echo "R=$R RUN_EXACT=$RUN_EXACT SKIP_EXISTING=$SKIP_EXISTING INCLUDE_GCOMP=$INCLUDE_GCOMP"
echo "Rout=$rout_file"

R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args ${r_args}" \
  "$run_file" \
  "$rout_file"
