#!/bin/bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/env.sh"
[[ "$#" -le 1 ]] || { echo 'Usage: bash submit.sh [comma-separated task IDs]'; exit 1; }
requested=${1:-all}
mkdir -p -- "$RRFULL_OUT/logs" "$RRFULL_OUT/.locks"
exec 9>"$RRFULL_OUT/.locks/submit.lock"
flock -n 9 || { echo 'Another submission is preparing this run'; exit 1; }
if [[ -f "$RRFULL_OUT/active_job_ids.txt" ]]; then
    while read -r id; do
        [[ "$id" =~ ^[0-9]+$ ]] || continue
        if [[ -n "$(squeue -h -j "$id" 2>/dev/null || true)" ]]; then
            echo "Formal RR job $id is still active; wait before resubmitting this run."; exit 1
        fi
    done < "$RRFULL_OUT/active_job_ids.txt"
fi
cd -- "$RRFULL_HOME"
sha256sum --quiet -c SHA256SUMS
module load StdEnv/2023 gcc/12.3 r/4.3.1
Rscript --vanilla prepare.R "$RRFULL_HOME" "$RRFULL_OUT" "$requested" | tee "$RRFULL_OUT/preflight.log"
ids=()
: > "$RRFULL_OUT/active_job_ids.txt"
while IFS=$'\t' read -r group cores mem walltime todo; do
    [[ -n "$group" ]] || continue
    [[ "$todo" =~ ^[0-9]+(,[0-9]+)*$ ]] || exit 1
    job=$(sbatch --parsable --account="$RRFULL_ACCOUNT" --job-name="RRv2_1000_g$group" \
        --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="$mem" --time="$walltime" \
        --array="$todo" --export=ALL --chdir="$RRFULL_HOME" \
        --output="$RRFULL_OUT/logs/RRfull_%A_%a.log" "$RRFULL_HOME/worker.sh")
    job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
    ids+=("$job");printf '%s\n' "$job" >> "$RRFULL_OUT/active_job_ids.txt"
    printf 'Submitted array %s: tasks=%s, %s CPUs, %s, %s; each task runs all 1000 seeds.\n' "$job" "$todo" "$cores" "$mem" "$walltime"
done < "$RRFULL_OUT/submission_groups.tsv"
if [[ ${#ids[@]} -eq 0 ]]; then echo 'All selected formal cells already processed; no simulation submitted.'; exit 0; fi
dependency=$(IFS=:; printf '%s' "${ids[*]}")
report=$(sbatch --parsable --account="$RRFULL_ACCOUNT" --job-name=RRv2_report --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=8G --time=01:00:00 --dependency="afterany:$dependency" --export=ALL \
    --chdir="$RRFULL_HOME" --output="$RRFULL_OUT/logs/report_%j.log" "$RRFULL_HOME/collect.sh")
report=${report%%;*}; [[ "$report" =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$report" >> "$RRFULL_OUT/active_job_ids.txt"
printf 'Report job: %s\nResults: %s\n' "$report" "$RRFULL_OUT"
