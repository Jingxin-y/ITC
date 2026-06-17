#!/bin/sh

export code_dir="/home/yanjin41/brmplus_simulation/compareITC/sandwich"
export result_dir="/scratch/yanjin41/RRRDOR/brmplus_simulation/ITC/sandwich"

module load StdEnv/2023
module load gcc/12.3 r/4.3.1

n=${1:?Usage: submit_simucpp_anchored.sh n R [event] [hypothesis] [run_exact] [ess_ratio]}
R=${2:?Usage: submit_simucpp_anchored.sh n R [event] [hypothesis] [run_exact] [ess_ratio]}
event=${3:-common}
hypothesis=${4:-alternative}
run_exact=${5:-TRUE}
ess_ratio=${6:-}


r_args="n=${n} R=${R} event='${event}' hypothesis='${hypothesis}' result_dir='${result_dir}' run_exact=${run_exact}"
out_suffix=""

if [ -n "$ess_ratio" ]; then
  r_args="${r_args} ess_ratio=${ess_ratio}"
  out_suffix="_essratio_${ess_ratio}"
fi

R --vanilla --max-connections=512 CMD BATCH --no-save --no-restore \
  "--args ${r_args}" \
  "$code_dir/run_simulation_anchored.R" \
  "$result_dir/Rout/RR_simucpp_anchored_${event}_${hypothesis}_N_${n}_R_${R}${out_suffix}.Rout"
