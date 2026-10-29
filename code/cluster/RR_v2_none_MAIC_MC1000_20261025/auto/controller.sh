#!/bin/bash
set -euo pipefail
round=${1:?Supply continuation round}
source "$RRNMW_AUTO_BASE/env.sh"
module load StdEnv/2023 gcc/12.3 r/4.3.1
cd -- "$RRNMW_AUTO_HOME"
sha256sum --quiet -c SHA256SUMS
(cd -- "$RRNMW_AUTO_BASE" && sha256sum --quiet -c SHA256SUMS)
exec 8>"$RRNMW_OUT/auto_none_MAIC_20261025/controller.lock"
flock -n 8 || { echo 'Another none controller is active'; exit 1; }
Rscript --vanilla "$RRNMW_AUTO_HOME/inspect_none.R" "$RRNMW_AUTO_BASE" "$RRNMW_OUT" "$RRNMW_AUTO_STATE" "$round" "$RRNMW_AUTO_MAX_ROUNDS"
action=$(cat "$RRNMW_AUTO_STATE/action.txt")
registry="$RRNMW_OUT/auto_none_MAIC_20261025/active_job_ids.txt"
remember_job() {
    [[ "$1" =~ ^[0-9]+$ ]] || exit 1
    printf '%s\n' "$1" >> "$registry"
    printf '%s\n' "$1" >> "$RRNMW_OUT/active_job_ids.txt"
}
if [[ "$action" == COMPLETE ]]; then
    echo 'All none cases have 1000 complete seed records. No further none computation will be queued.'
    # Refresh the full report after any remaining original score/report jobs.
    remaining=()
    queue=$(squeue -h --array -u "$USER" -o '%F|%j')
    while read -r id; do
        [[ "$id" =~ ^[0-9]+$ && "$id" != "${SLURM_JOB_ID:-}" ]] || continue
        case " $(cat "$registry" | tr '\n' ' ') " in *" $id "*) continue ;; esac
        if awk -F'|' -v id="$id" '$1==id {found=1} END{exit !found}' <<< "$queue"; then remaining+=("$id");fi
    done < "$RRNMW_OUT/active_job_ids.txt"
    dep=()
    if [[ ${#remaining[@]} -gt 0 ]]; then joined=$(IFS=:; printf '%s' "${remaining[*]}");dep=(--dependency="afterany:$joined");fi
    report=$(sbatch --parsable --account="$RRNMW_ACCOUNT" --job-name=RRnoneMAIC_finalreport --nodes=1 --ntasks=1 \
        --cpus-per-task=1 --mem=8G --time=01:00:00 "${dep[@]}" --export=ALL --chdir="$RRNMW_AUTO_BASE" \
        --output="$RRNMW_AUTO_STATE/final_report_%j.log" "$RRNMW_AUTO_BASE/collect.sh")
    report=${report%%;*};remember_job "$report"
    printf 'COMPLETE\n' > "$RRNMW_AUTO_STATE/status.txt";printf '%s\n' "$report" > "$RRNMW_AUTO_STATE/final_report_job_id.txt"
    echo "Final full-report refresh queued: $report";exit 0
fi
if [[ "$action" != RESUME ]]; then
    printf '%s\n' "$action" > "$RRNMW_AUTO_STATE/status.txt"
    echo "Automatic none continuation stopped: $action. See pending_seeds.csv and latest_summary.csv."
    exit 2
fi
todo=$(cat "$RRNMW_AUTO_STATE/pending_tasks.txt")
[[ "$todo" =~ ^[0-9]+(,[0-9]+)*$ ]] || exit 1
IFS=, read -ra ids <<< "$todo"
for id in "${ids[@]}";do case "$id" in 1|2|3|4|5|6|7|8) ;; *) echo 'Non-none task rejected';exit 1;;esac;done
job=$(sbatch --parsable --account="$RRNMW_ACCOUNT" --job-name=RRnoneMAIC_resume --nodes=1 --ntasks=1 \
    --cpus-per-task=64 --mem=10G --time=04:00:00 --array="$todo" --export=ALL --chdir="$RRNMW_AUTO_BASE" \
    --output="$RRNMW_AUTO_STATE/resume_%A_%a.log" "$RRNMW_AUTO_BASE/worker.sh")
job=${job%%;*};remember_job "$job"
next=$((round+1))
check=$(sbatch --parsable --account="$RRNMW_ACCOUNT" --job-name=RRnoneMAIC_check --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=4G --time=00:30:00 --dependency="afterany:$job" --export=ALL --chdir="$RRNMW_AUTO_HOME" \
    --output="$RRNMW_AUTO_STATE/check_%j.log" "$RRNMW_AUTO_HOME/controller.sh" "$next")
check=${check%%;*};remember_job "$check"
printf 'RESUME_ROUND_%s\n' "$next" > "$RRNMW_AUTO_STATE/status.txt"
echo "None continuation round $next queued: array $job tasks $todo; next checker $check"
