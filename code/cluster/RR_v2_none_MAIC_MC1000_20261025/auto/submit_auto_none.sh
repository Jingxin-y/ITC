#!/bin/bash
set -euo pipefail
auto_home=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export RRNMW_AUTO_HOME="$auto_home"
export RRNMW_AUTO_BASE="${RRNMW_AUTO_BASE:-$(cd -- "$auto_home/.." && pwd -P)}"
export RRNMW_HOME="$RRNMW_AUTO_BASE"
source "$RRNMW_AUTO_BASE/env.sh"
export RRNMW_AUTO_MAX_ROUNDS="${RRNMW_AUTO_MAX_ROUNDS:-10}"
[[ "$#" -le 1 && "$RRNMW_AUTO_MAX_ROUNDS" =~ ^[1-9][0-9]*$ ]] || { echo 'Usage: bash submit_auto_none.sh [current_none_array_job_id]'; exit 1; }
[[ -f "$RRNMW_OUT/active_job_ids.txt" && -f "$RRNMW_OUT/submission_groups.tsv" ]] || { echo "Original run metadata not found: $RRNMW_OUT"; exit 1; }
mkdir -p -- "$RRNMW_OUT/auto_none_MAIC_20261025" "$RRNMW_OUT/logs" "$RRNMW_OUT/.locks"
exec 9>"$RRNMW_OUT/.locks/submit.lock"
flock -n 9 || { echo 'Another submission is preparing this result directory'; exit 1; }
registry="$RRNMW_OUT/auto_none_MAIC_20261025/active_job_ids.txt"
queue=$(squeue -h --array -u "$USER" -o '%F|%j')
if [[ -f "$registry" ]]; then
    while read -r id; do
        [[ "$id" =~ ^[0-9]+$ ]] || continue
        if awk -F'|' -v id="$id" '$1==id {found=1} END{exit !found}' <<< "$queue"; then
            echo "None automatic continuation is already queued/running: $id"; exit 0
        fi
    done < "$registry"
fi
cd -- "$auto_home"
sha256sum --quiet -c SHA256SUMS
# Query the queue once. Match both this run's recorded IDs and groups consisting
# solely of the none task IDs; score jobs and other result directories are excluded.
parents=()
if [[ "$#" == 1 ]]; then
    [[ "$1" =~ ^[0-9]+$ ]] || exit 1
    if awk -F'|' -v id="$1" '$1==id {found=1} END{exit !found}' <<< "$queue"; then parents+=("$1"); fi
else
    while IFS=$'\t' read -r group cores memory walltime ids; do
        [[ -n "$group" ]] || continue
        is_none=1
        IFS=, read -ra values <<< "$ids"
        for task in "${values[@]}"; do
            case "$task" in 1|2|3|4|5|6|7|8) ;; *) is_none=0 ;; esac
        done
        [[ "$is_none" == 1 ]] || continue
        while IFS='|' read -r id name; do
            [[ "$name" == "RRnoneMAIC_g$group" && "$id" =~ ^[0-9]+$ ]] || continue
            grep -qx "$id" "$RRNMW_OUT/active_job_ids.txt" || continue
            parents+=("$id")
        done <<< "$queue"
    done < "$RRNMW_OUT/submission_groups.tsv"
fi
if [[ ${#parents[@]} -gt 0 ]]; then mapfile -t parents < <(printf '%s\n' "${parents[@]}" | sort -u); fi
export RRNMW_AUTO_STATE="$RRNMW_OUT/auto_none_MAIC_20261025/chain_$(date +%Y%m%d_%H%M%S)_$$"
mkdir -p -- "$RRNMW_AUTO_STATE"
printf '%s\n' "$RRNMW_AUTO_STATE" > "$RRNMW_OUT/auto_none_MAIC_20261025/latest_chain.txt"
dependency=()
if [[ ${#parents[@]} -gt 0 ]]; then
    joined=$(IFS=:; printf '%s' "${parents[*]}");dependency=(--dependency="afterany:$joined")
    echo "Check will wait for the current none array(s): $joined"
else echo 'No active none array found; queue the check now. It will verify all 8000 completion flags.'; fi
job=$(sbatch --parsable --account="$RRNMW_ACCOUNT" --job-name=RRnoneMAIC_check --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=4G --time=00:30:00 "${dependency[@]}" --export=ALL --chdir="$auto_home" \
    --output="$RRNMW_AUTO_STATE/check_%j.log" "$auto_home/controller.sh" 0)
job=${job%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || exit 1
printf '%s\n' "$job" >> "$registry"
printf '%s\n' "$job" >> "$RRNMW_OUT/active_job_ids.txt"
printf 'Queued none checker: %s\nState and progress: %s\n' "$job" "$RRNMW_AUTO_STATE"
