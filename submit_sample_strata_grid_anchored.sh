#!/bin/sh
#SBATCH --job-name=strata_grid
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=8G

set -eu

export code_dir="${code_dir:-/home/yanjin41/brmplus_simulation/compareITC/sandwich}"
export result_dir="${result_dir:-/scratch/yanjin41/RRRDOR/brmplus_simulation/ITC/sandwich}"
export output_dir="${output_dir:-${result_dir}/sample_strata_grid}"

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

mkdir -p "$output_dir" "$result_dir/Rout"

run_file="$code_dir/summarize_sample_strata_grid_anchored.R"
job_id="${SLURM_JOB_ID:-manual}"
rout_file="$result_dir/Rout/RR_sample_strata_anchored_grid_R_1000_${job_id}.Rout"

if [ ! -f "$run_file" ]; then
  echo "Cannot find run file: $run_file" >&2
  exit 1
fi

r_args="R=1000 first_r=1 last_r=1000 result_dir='${result_dir}' output_dir='${output_dir}' code_dir='${code_dir}'"

echo "Summarizing anchored sample strata grid"
echo "run_file=$run_file"
echo "result_dir=$result_dir"
echo "output_dir=$output_dir"
echo "Rout=$rout_file"

R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args ${r_args}" \
  "$run_file" \
  "$rout_file"
