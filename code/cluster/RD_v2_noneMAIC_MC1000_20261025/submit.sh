#!/bin/bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/env.sh"
[[ "$#" -le 1 ]] || { echo 'Usage: bash submit.sh [comma-separated task IDs]'; exit 1; }
requested=${1:-all}
mkdir -p -- "$RDNM_OUT/logs" "$RDNM_OUT/.locks"
exec 9>"$RDNM_OUT/.locks/submit.lock"
flock -n 9 || { echo 'Another submission is preparing this run'; exit 1; }
if [[ -f "$RDNM_OUT/active_job_ids.txt" ]]; then
    while read -r id; do
        [[ "$id" =~ ^[0-9]+$ ]] || continue
        if [[ -n "$(squeue -h -j "$id" 2>/dev/null || true)" ]]; then
            echo "Formal RD job $id is still active; wait before resubmitting this run."; exit 1
        fi
    done < "$RDNM_OUT/active_job_ids.txt"
fi
cd -- "$RDNM_HOME"
sha256sum --quiet -c SHA256SUMS
module load StdEnv/2023 gcc/12.3 r/4.3.1
Rscript --vanilla prepare.R "$RDNM_HOME" "$RDNM_OUT" "$requested" | tee "$RDNM_OUT/preflight.log"
ids=()
: > "$RDNM_OUT/active_job_ids.txt"
while IFS=$'\t' read -r group cores mem walltime todo; do
    [[ -n "$group" ]] || continue
    [[ "$todo" =~ ^[0-9]+(,[0-9]+)*$ ]] || exit 1
    job=$(sbatch --parsable --account="$RDNM_ACCOUNT" --job-name="RDnoneMAIC_g$group" \
        --nodes=1 --ntasks=1 --cpus-per-task="$cores" --mem="$mem" --time="$walltime" \
        --array="$todo" --export=ALL --chdir="$RDNM_HOME" \
        --output="$RDNM_OUT/logs/RDnoneMAIC_%A_%a.log" "$RDNM_HOME/worker.sh")
    job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
    ids+=("$job");printf '%s\n' "$job" >> "$RDNM_OUT/active_job_ids.txt"
    printf 'Submitted array %s: tasks=%s, %s CPUs, %s, %s; each task runs all 1000 seeds.\n' "$job" "$todo" "$cores" "$mem" "$walltime"
done < "$RDNM_OUT/submission_groups.tsv"
if [[ ${#ids[@]} -eq 0 ]]; then echo 'All selected formal cells already processed; no simulation submitted.'; exit 0; fi
dependency=$(IFS=:; printf '%s' "${ids[*]}")
report=$(sbatch --parsable --account="$RDNM_ACCOUNT" --job-name=RDnoneMAIC_report --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=8G --time=01:00:00 --dependency="afterany:$dependency" --export=ALL \
    --chdir="$RDNM_HOME" --output="$RDNM_OUT/logs/report_%j.log" "$RDNM_HOME/collect.sh")
report=${report%%;*}; [[ "$report" =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$report" >> "$RDNM_OUT/active_job_ids.txt"
printf 'Report job: %s\nResults: %s\n' "$report" "$RDNM_OUT"
