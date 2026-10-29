#!/bin/bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/env.sh"
[[ "$#" -le 1 ]] || { echo 'Usage: bash submit.sh [comma-separated task IDs]'; exit 1; }
requested=${1:-all}
mkdir -p -- "$RDFULL_OUT/logs" "$RDFULL_OUT/.locks"
exec 9>"$RDFULL_OUT/.locks/submit.lock"
flock -n 9 || { echo 'Another submission is preparing this run'; exit 1; }
if [[ -f "$RDFULL_OUT/active_job_ids.txt" ]]; then
    while read -r id; do
        [[ "$id" =~ ^[0-9]+$ ]] || continue
        if [[ -n "$(squeue -h -j "$id" 2>/dev/null || true)" ]]; then
            echo "Formal RD job $id is still active; wait before resubmitting this run."; exit 1
        fi
    done < "$RDFULL_OUT/active_job_ids.txt"
fi
cd -- "$RDFULL_HOME"
sha256sum --quiet -c SHA256SUMS
module load StdEnv/2023 gcc/12.3 r/4.3.1
Rscript --vanilla prepare.R "$RDFULL_HOME" "$RDFULL_OUT" "$requested" | tee "$RDFULL_OUT/preflight.log"
ids=()
: > "$RDFULL_OUT/active_job_ids.txt"
while IFS=$'\t' read -r group cores mem walltime todo; do
    [[ -n "$group" ]] || continue
    [[ "$todo" =~ ^[0-9]+(,[0-9]+)*$ ]] || exit 1
    job=$(sbatch --parsable --account="$RDFULL_ACCOUNT" --job-name="RDv2_1000_g$group" \
        --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="$mem" --time="$walltime" \
        --array="$todo" --export=ALL --chdir="$RDFULL_HOME" \
        --output="$RDFULL_OUT/logs/RDfull_%A_%a.log" "$RDFULL_HOME/worker.sh")
    job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
    ids+=("$job");printf '%s\n' "$job" >> "$RDFULL_OUT/active_job_ids.txt"
    printf 'Submitted array %s: tasks=%s, %s CPUs, %s, %s; each task runs all 1000 seeds.\n' "$job" "$todo" "$cores" "$mem" "$walltime"
done < "$RDFULL_OUT/submission_groups.tsv"
if [[ ${#ids[@]} -eq 0 ]]; then echo 'All selected formal cells already processed; no simulation submitted.'; exit 0; fi
dependency=$(IFS=:; printf '%s' "${ids[*]}")
report=$(sbatch --parsable --account="$RDFULL_ACCOUNT" --job-name=RDv2_report --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=8G --time=01:00:00 --dependency="afterany:$dependency" --export=ALL \
    --chdir="$RDFULL_HOME" --output="$RDFULL_OUT/logs/report_%j.log" "$RDFULL_HOME/collect.sh")
report=${report%%;*}; [[ "$report" =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$report" >> "$RDFULL_OUT/active_job_ids.txt"
printf 'Report job: %s\nResults: %s\n' "$report" "$RDFULL_OUT"
