#!/usr/bin/env bash
set -uo pipefail

repo_root="${1:-/home/mrk/sbst/evo-subster}"
cd "$repo_root" || { echo "cannot enter repository: $repo_root" >&2; exit 2; }

# Load the same noninteractive tool modules used by the interactive biohazard shell.
if [[ -f /etc/profile.d/lmod.sh ]]; then
    source /etc/profile.d/lmod.sh
    module use /big/mrk/app/.modulefiles
    module load last/1648 yq/4.45.1 jq/1.8.1 datasets/18.23.0 dataformat/18.9.0
fi
for required_tool in lastal lastdb last-train yq jq Rscript; do
    command -v "$required_tool" >/dev/null 2>&1 || { echo "missing required tool: $required_tool" >&2; exit 2; }
done

date_tag="20260818"
result_root="results/glires_20260818"
selection_tsv="$result_root/trio_selection/selected_trios.tsv"
cache_dir="$result_root/trio_selection/train_cache"
log_dir="$result_root/logs"
run_id="glires_no_busco_94g_p1_$(date +%Y%m%dT%H%M%S)"
summary="$log_dir/${run_id}_summary.tsv"
status_file="$log_dir/${run_id}.exit_status"

export LASTAL_QUERY_BATCH_SIZE=8K
export LASTAL_SPLIT_MEMORY=94G
export THREAD_NUM=1

mkdir -p "$log_dir"
[[ -s "$selection_tsv" && -d "$cache_dir" ]] || {
    echo "missing selected-trio table or train cache" >&2
    exit 2
}
if command -v flock >/dev/null 2>&1; then
    exec 9>"$log_dir/.glires_no_busco_94g_p1.lock"
    flock -n 9 || { echo "another Glires runner already holds the lock" >&2; exit 3; }
fi
printf 'trio\tout_acc\tin1_acc\tin2_acc\tstatus\tattempts\tlog\tnote\n' >"$summary"

tsv_field() {
    local header="$1" name="$2"
    awk -F '\t' -v target="$name" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == target) { print i; exit } }
    ' "$header"
}

out_col=$(tsv_field "$selection_tsv" out_acc)
in1_col=$(tsv_field "$selection_tsv" in1_acc)
in2_col=$(tsv_field "$selection_tsv" in2_acc)
out_short_col=$(tsv_field "$selection_tsv" out_short)
in1_short_col=$(tsv_field "$selection_tsv" in1_short)
in2_short_col=$(tsv_field "$selection_tsv" in2_short)
for col in "$out_col" "$in1_col" "$in2_col" "$out_short_col" "$in1_short_col" "$in2_short_col"; do
    [[ "$col" =~ ^[0-9]+$ ]] || { echo "selected_trios.tsv is missing required columns" >&2; exit 2; }
done

mem_available_gib() {
    awk '/MemAvailable:/ { printf "%d", $2 / 1048576 }' /proc/meminfo
}

load_1m() {
    awk '{ print $1 }' /proc/loadavg
}

swap_used_gib() {
    free -g | awk '/^Swap:/ { print $3 }'
}

wait_for_shared_resources() {
    local available load swap_used
    while true; do
        available=$(mem_available_gib)
        load=$(load_1m)
        swap_used=$(swap_used_gib)
        if (( available >= 200 && swap_used <= 4 )) && awk -v l="$load" 'BEGIN { exit !(l < 32) }'; then
            return 0
        fi
        printf '%s waiting: MemAvailable=%sGiB SwapUsed=%sGiB load1=%s\n' "$(date -Is)" "$available" "$swap_used" "$load" >&2
        sleep 60 9>&-
    done
}

is_complete_run() {
    local out_acc="$1" in1_acc="$2" in2_acc="$3"
    python3 - "$result_root" "$date_tag" "$out_acc" "$in1_acc" "$in2_acc" <<'PYRUN'
import sys
from pathlib import Path
sys.path.insert(0,str(Path('src/lib').resolve()))
import run_manifest as runs
root=Path(sys.argv[1])
if not root.exists(): raise SystemExit(1)
try:
    found=runs.discover(root)
except ValueError:
    raise SystemExit(1)
for run in found:
    manifest=runs.load(run)
    if manifest['date']==sys.argv[2] and [x['accession'] for x in manifest['organisms']]==sys.argv[3:]:
        runs.validated(run,complete=True)
        raise SystemExit(0)
raise SystemExit(1)
PYRUN
}

is_memory_failure() {
    grep -Eqi 'skipping sequence|LASTAL_SPLIT_MEMORY=.*too low|Cannot allocate memory|out of memory|oom-kill' "$1"
}

run_trio() {
    local trio="$1" out_acc="$2" in1_acc="$3" in2_acc="$4"
    local out_short="$5" in1_short="$6" in2_short="$7"
    local dir log rc attempts=0 note
    log="$log_dir/${run_id}_${trio}.log"


    while (( attempts < 2 )); do
        attempts=$((attempts + 1))
        wait_for_shared_resources
        printf '%s attempt=%s trio=%s MemAvailable=%sGiB load1=%s\n' "$(date -Is)" "$attempts" "$trio" "$(mem_available_gib)" "$(load_1m)" >>"$log"
        if command -v ionice >/dev/null 2>&1; then
            nice -n 10 ionice -c 2 -n 7 bash src/sbst_fromDwl.sh "$date_tag" "$out_acc" "$in1_acc" "$in2_acc" \
                --genome-dir genomes --out-dir "$result_root" --train-cache-dir "$cache_dir" --thread 1 >>"$log" 2>&1
        else
            nice -n 10 bash src/sbst_fromDwl.sh "$date_tag" "$out_acc" "$in1_acc" "$in2_acc" \
                --genome-dir genomes --out-dir "$result_root" --train-cache-dir "$cache_dir" --thread 1 >>"$log" 2>&1
        fi
        rc=$?
        if is_memory_failure "$log"; then
            printf '%s\t%s\t%s\t%s\tskipped_memory_cap\t%s\t%s\tLAST memory-cap diagnostic\n' "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$attempts" "$log" >>"$summary"
            return 0
        fi
        if [[ "$rc" -eq 0 ]] && is_complete_run "$out_acc" "$in1_acc" "$in2_acc"; then
            printf '%s\t%s\t%s\t%s\tcompleted\t%s\t%s\tvalidated\n' "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$attempts" "$log" >>"$summary"
            return 0
        fi
        if (( attempts == 1 )); then
            printf '%s retrying through manifest-verified stage recovery\n' "$(date -Is)" >>"$log"
        fi
    done
    note="pipeline exit=${rc:-unknown} or final-artifact validation failed"
    printf '%s\t%s\t%s\t%s\tfailed_other\t%s\t%s\t%s\n' "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$attempts" "$log" "$note" >>"$summary"
}

while IFS=$'\t' read -r out_acc in1_acc in2_acc out_short in1_short in2_short; do
    [[ -n "$out_acc" && -n "$in1_acc" && -n "$in2_acc" ]] || continue
    trio="${out_short}1_${in1_short}2_${in2_short}3"
    run_trio "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$out_short" "$in1_short" "$in2_short"
done < <(awk -F '\t' -v o="$out_col" -v a="$in1_col" -v b="$in2_col" -v os="$out_short_col" -v as="$in1_short_col" -v bs="$in2_short_col" 'NR > 1 { print $o "\t" $a "\t" $b "\t" $os "\t" $as "\t" $bs }' "$selection_tsv")

printf '0\n' >"$status_file"
printf 'summary=%s\n' "$summary"
