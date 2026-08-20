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

run_dir_for() {
    printf '%s/%s1_%s2_%s3/%s\n' "$result_root" "$1" "$2" "$3" "$date_tag"
}

is_complete_run() {
    local dir="$1" out_short="$2" in1_short="$3" in2_short="$4"
    local out_acc="$5" in1_acc="$6" in2_acc="$7"
    local out_label="${out_short}1" in1_label="${in1_short}2" in2_label="${in2_short}3"
    local joined="$dir/intermediateFiles/${out_label}_${in1_label}_${in2_label}_${date_tag}.maf"
    local ncds="${joined%.maf}_ncds.maf"
    local ratio="$dir/statistics/misc/sbstRatio_${date_tag}.out"
    local short acc kind suffix stat pdf manifest requires_ncds=0

    [[ -s "$joined" && -e "${joined}.complete" && -s "$ratio" && -e "${ratio}.complete" ]] || return 1
    manifest="$dir/metadata/metadata_manifest.json"
    if [[ -s "$manifest" ]] && jq -e 'any(.organisms[]; .slot == "org1" and (.gff_path != null and .gff_path != ""))' "$manifest" >/dev/null; then
        requires_ncds=1
        [[ -s "$ncds" ]] || return 1
    fi
    for short in "$in1_label" "$in2_label"; do
        if [[ "$short" == "$in1_label" ]]; then acc="$in1_acc"; else acc="$in2_acc"; fi
        for kind in singlenuc dinuc; do
            for suffix in '' $([[ "$requires_ncds" -eq 1 ]] && printf '_ncds'); do
                stat="$dir/statistics/$short/$kind/${acc}_${short}_${date_tag}"
                if [[ "$kind" == "dinuc" ]]; then stat+="_dinuc"; fi
                stat+="${suffix}.tsv"
                [[ -s "$stat" ]] || return 1
            done
        done
    done
    for kind in singlenuc dinuc; do
        for suffix in '' $([[ "$requires_ncds" -eq 1 ]] && printf '_ncds'); do
            stat="$dir/statistics/$in1_label/$kind/${in1_acc}_${in1_label}_${date_tag}"
            [[ "$kind" == "dinuc" ]] && stat+="_dinuc"
            stat+="${suffix}.tsv.complete"
            [[ -e "$stat" ]] || return 1
        done
    done
    for short in "$in1_label" "$in2_label"; do
        if [[ "$short" == "$in1_label" ]]; then acc="$in1_acc"; else acc="$in2_acc"; fi
        for pdf in \
            "figs/$short/dinuc/${acc}_${short}_${date_tag}_dinuc.tsv.pdf" \
            "figs/$short/singlenuc/count/${acc}_${short}_${date_tag}_ori.pdf" \
            "figs/$short/singlenuc/count/${acc}_${short}_${date_tag}_sbst.pdf" \
            "figs/$short/singlenuc/log-ratio/${acc}_${short}_${date_tag}_logRatio.pdf" \
            "figs/$short/singlenuc/ratio/${acc}_${short}_${date_tag}_norm.pdf"; do
            [[ -s "$dir/$pdf" ]] || return 1
        done
        if [[ "$requires_ncds" -eq 1 ]]; then
            for pdf in \
                "figs/$short/dinuc/${acc}_${short}_${date_tag}_dinuc_ncds.tsv.pdf" \
                "figs/$short/singlenuc/count/${acc}_${short}_${date_tag}_ncds_ori.pdf" \
                "figs/$short/singlenuc/count/${acc}_${short}_${date_tag}_ncds_sbst.pdf" \
                "figs/$short/singlenuc/log-ratio/${acc}_${short}_${date_tag}_ncds_logRatio.pdf" \
                "figs/$short/singlenuc/ratio/${acc}_${short}_${date_tag}_ncds_norm.pdf"; do
                [[ -s "$dir/$pdf" ]] || return 1
            done
        fi
    done
}

safe_remove_run_dir() {
    local dir="$1" expected="$2" base_real dir_real
    base_real=$(realpath -m "$result_root") || return 1
    dir_real=$(realpath -m "$dir") || return 1
    [[ "$dir_real" == "$base_real/"* && "$dir_real" == "$expected" ]] || {
        echo "refusing unsafe removal: $dir_real" >&2
        return 1
    }
    [[ ! -L "$dir" ]] || { echo "refusing symlink run dir: $dir" >&2; return 1; }
    rm -rf -- "$dir"
}

is_memory_failure() {
    grep -Eqi 'skipping sequence|LASTAL_SPLIT_MEMORY=.*too low|Cannot allocate memory|out of memory|oom-kill' "$1"
}

run_trio() {
    local trio="$1" out_acc="$2" in1_acc="$3" in2_acc="$4"
    local out_short="$5" in1_short="$6" in2_short="$7"
    local dir log rc attempts=0 note
    dir=$(run_dir_for "$out_short" "$in1_short" "$in2_short")
    log="$log_dir/${run_id}_${trio}.log"

    if is_complete_run "$dir" "$out_short" "$in1_short" "$in2_short" "$out_acc" "$in1_acc" "$in2_acc"; then
        printf '%s\t%s\t%s\t%s\texisting_complete\t0\t%s\tvalidated\n' "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$log" >>"$summary"
        return 0
    fi

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
        if [[ "$rc" -eq 0 ]] && is_complete_run "$dir" "$out_short" "$in1_short" "$in2_short" "$out_acc" "$in1_acc" "$in2_acc"; then
            printf '%s\t%s\t%s\t%s\tcompleted\t%s\t%s\tvalidated\n' "$trio" "$out_acc" "$in1_acc" "$in2_acc" "$attempts" "$log" >>"$summary"
            return 0
        fi
        if (( attempts == 1 )); then
            expected="$(realpath -m "$result_root")/${out_short}1_${in1_short}2_${in2_short}3/${date_tag}"
            printf '%s resume validation failed; removing only this incomplete run before full retry\n' "$(date -Is)" >>"$log"
            safe_remove_run_dir "$dir" "$expected" || break
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
