#!/bin/bash

# Shared, shell-free LAST option parsing (also inherited by tree-mode children).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$ROOT_DIR/src/lib/last_options.sh"
last_options_parse_cli "$@" || exit 1
set -- "${LAST_REMAINING_ARGS[@]}"

for _arg in "$@"; do
    if [[ "$_arg" == "-h" || "$_arg" == "--help" ]]; then
        cat <<'EOF'
Usage: sbst.sh <DATE> <ORG1.fa> <ORG2.fa> <ORG3.fa> <ORG1.gff|NO_GFF_FILE> [OPTIONS]

Run the substitution spectrum pipeline from existing FASTA files.

Positional arguments:
  DATE                    Run date (e.g. 20240101)
  ORG1.fa                 Outgroup FASTA file path
  ORG2.fa                 Ingroup FASTA file path
  ORG3.fa                 Ingroup FASTA file path
  ORG1.gff|NO_GFF_FILE    GFF annotation for outgroup, or literal NO_GFF_FILE

Options:
  --out-dir PATH              Output directory (default: ./results)
  --org2-gff PATH             Ingroup 2 GFF (default: NO_GFF_FILE)
  --org3-gff PATH             Ingroup 3 GFF (default: NO_GFF_FILE)
  --train-cache-dir PATH      Directory containing cached last-train files
  --accession-ids A1 A2 A3   Exact NCBI accessions for cache keys; normally supplied by sbst_fromDwl.sh
  --thread N                  Number of threads for LAST alignment (default: 8)
  --lastdb-options TEXT    Additional lastdb arguments (training and alignment)
  --lastal-options TEXT    Additional lastal arguments
  --last-train-options TEXT  Additional last-train arguments (including tree selection)
  --last-split-options TEXT  Additional last-split arguments
  --idt-only                  Stop after checking sequence percent identity among three genomes; skip downstream analysis
  -h, --help                  Show this help message and exit

LAST option strings use shell-style quoting without shell expansion (LAST 1648).
Specify each --*-options once. Pipeline-managed formats, inputs and threads
cannot be replaced. Different settings require a new date/output/cache directory.
EOF
        exit 0
    fi
done

lastal --version

# Resolve script locations and configuration relative to this script
ALIGN_DIR="$ROOT_DIR/src/align"
COUNT_DIR="$ROOT_DIR/src/count"
METRICS_DIR="$ROOT_DIR/src/metrics"
R_DIR="$ROOT_DIR/src/visualize"

# Allow callers to override directories, otherwise rely on inferred paths
ALIGN_DIR="${ALIGN_DIR_OVERRIDE:-$ALIGN_DIR}"
COUNT_DIR="${COUNT_DIR_OVERRIDE:-$COUNT_DIR}"
METRICS_DIR="${METRICS_DIR_OVERRIDE:-$METRICS_DIR}"
R_DIR="${R_DIR_OVERRIDE:-$R_DIR}"

# Ensure helper scripts are discoverable without absolute paths
PATH="$SCRIPT_DIR:$ALIGN_DIR:$PATH"

# shellcheck source=lib/train_cache.sh
source "$SCRIPT_DIR/lib/train_cache.sh"

config_file="$ROOT_DIR/config/sbst_config.yaml"
dwl_config_file="$ROOT_DIR/config/dwl_config.yaml"

# Load YAML configuration using yq
if [ ! -f "$config_file" ]; then
    echo "Configuration file not found!" 1>&2
    exit 1
fi
if [ ! -f "$dwl_config_file" ]; then
    echo "Download configuration file not found!" 1>&2
    exit 1
fi

# Function to get config values using yq
get_config() {
    yq eval "$1" "$config_file"
}

get_dwl_config() {
    yq eval "$1" "$dwl_config_file"
}

# Resolve a config path relative to the evo-subster repo root (ROOT_DIR)
# Leaves absolute paths untouched.
resolve_path() {
    local p="$1"
    if [[ "$p" != /* ]]; then
        p="${p#./}"
        p="$ROOT_DIR/$p"
    fi
    echo "$p"
}

OUT_DIR_OVERRIDE=""
TRAIN_CACHE_DIR=""
ACCESSION_IDS=()
IDT_ONLY=0
THREAD_NUM_OVERRIDE=8
org2GFF=""
org3GFF=""
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --org2-gff|--org3-gff|--org2-gff=*|--org3-gff=*)
            gff_option="${1%%=*}"
            if [[ "$1" == *=* ]]; then
                gff_value="${1#*=}"
                shift
            else
                gff_value="${2:-}"
                shift
                if [ "$#" -gt 0 ]; then shift; fi
            fi
            if [[ -z "$gff_value" || "$gff_value" == -* ]]; then
                echo "Error: $gff_option requires a non-empty GFF path." >&2
                exit 1
            fi
            if [ "$gff_option" = "--org2-gff" ]; then
                org2GFF="$gff_value"
            else
                org3GFF="$gff_value"
            fi
            continue
            ;;
        --idt-only)
            IDT_ONLY=1
            shift
            continue
            ;;
        --thread)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --thread requires a non-empty integer argument." >&2
                exit 1
            fi
            if ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]]; then
                echo "Error: --thread must be a positive integer (got: $2)." >&2
                exit 1
            fi
            THREAD_NUM_OVERRIDE="$2"
            shift 2
            continue
            ;;
        --thread=*)
            THREAD_NUM_OVERRIDE="${1#*=}"
            if [[ -z "$THREAD_NUM_OVERRIDE" ]]; then
                echo "Error: --thread requires a non-empty integer argument." >&2
                exit 1
            fi
            if ! [[ "$THREAD_NUM_OVERRIDE" =~ ^[0-9]+$ ]] || [[ "$THREAD_NUM_OVERRIDE" -lt 1 ]]; then
                echo "Error: --thread must be a positive integer (got: $THREAD_NUM_OVERRIDE)." >&2
                exit 1
            fi
            shift
            continue
            ;;
        --out-dir)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --out-dir requires a non-empty path argument." >&2
                exit 1
            fi
            OUT_DIR_OVERRIDE="$2"
            shift 2
            continue
            ;;
        --out-dir=*)
            OUT_DIR_OVERRIDE="${1#*=}"
            if [[ -z "$OUT_DIR_OVERRIDE" ]]; then
                echo "Error: --out-dir requires a non-empty path argument." >&2
                exit 1
            fi
            shift
            continue
            ;;
        --train-cache-dir)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --train-cache-dir requires a non-empty path argument." >&2
                exit 1
            fi
            TRAIN_CACHE_DIR="$2"
            shift 2
            continue
            ;;
        --train-cache-dir=*)
            TRAIN_CACHE_DIR="${1#*=}"
            if [[ -z "$TRAIN_CACHE_DIR" ]]; then
                echo "Error: --train-cache-dir requires a non-empty path argument." >&2
                exit 1
            fi
            shift
            continue
            ;;
        --accession-ids)
            if [ "$#" -lt 4 ]; then
                echo "Error: --accession-ids requires exactly three NCBI accessions." >&2
                exit 1
            fi
            for accession in "$2" "$3" "$4"; do
                if ! [[ "$accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
                    echo "Error: --accession-ids values must be versioned GCA/GCF accessions (got: $accession)." >&2
                    exit 1
                fi
            done
            ACCESSION_IDS=("$2" "$3" "$4")
            shift 4
            continue
            ;;
        --)
            shift
            POSITIONAL_ARGS+=("$@")
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL_ARGS[@]}"

if [ -n "$TRAIN_CACHE_DIR" ]; then
    TRAIN_CACHE_DIR=$(resolve_path "$TRAIN_CACHE_DIR")
fi

extract_accession_from_path() {
    local input_path=$1
    local base
    base=$(basename "$input_path")

    if [[ -z $base ]]; then
        echo ""
        return 0
    fi

    local first_part="${base%%_*}"
    local after_first="${base#*_}"

    if [[ $after_first == "$base" ]]; then
        echo "$base"
        return 0
    fi

    local second_part="${after_first%%_*}"

    if [[ -z $second_part ]]; then
        echo "$first_part"
        return 0
    fi

    echo "${first_part}_${second_part}"
}

make_short_name() {
    local full_name="$1"
    local suffix="$2"
    local first_part=""
    local second_part=""
    local rest=""
    local IFS='_'
    read -r first_part second_part rest <<< "$full_name"

    local first_trimmed="${first_part:0:3}"

    if [ -z "$second_part" ]; then
        second_part="${first_part:3}"
    fi

    local second_trimmed="${second_part:0:3}"

    echo "${first_trimmed}${second_trimmed}${suffix}"
}

# Get required arguments count from config
argNum=$(get_config '.settings.required_args')
if [ $# -ne "$argNum" ]; then
    echo "$(get_config '.errors.arg_count' | sed "s/{arg_num}/$argNum/g")" 1>&2
    echo "$(get_config '.errors.usage')" 1>&2
    exit 1
fi

DATE=$1
org1FASTA=$2
org2FASTA=$3
org3FASTA=$4
org1GFF=$5

# Resolve annotations before changing into the output directory. The outgroup
# sentinel remains authoritative; only omitted ingroup paths are auto-detected.
gff_fastas=("$org1FASTA" "$org2FASTA" "$org3FASTA")
gff_paths=("$org1GFF" "$org2GFF" "$org3GFF")
for i in 0 1 2; do
    if [ "$i" -gt 0 ] && [ -z "${gff_paths[$i]}" ]; then
        gff_paths[$i]="NO_GFF_FILE"
    fi
    if [ "${gff_paths[$i]}" != "NO_GFF_FILE" ]; then
        if [ ! -f "${gff_paths[$i]}" ] || [ ! -r "${gff_paths[$i]}" ]; then
            echo "Error: GFF file is not readable: ${gff_paths[$i]}" >&2
            exit 1
        fi
        gff_paths[$i]="$(cd "$(dirname "${gff_paths[$i]}")" && pwd)/$(basename "${gff_paths[$i]}")"
    fi
done
org1GFF="${gff_paths[0]}"

# Propagate thread count to child scripts (last_train.sh / one2one.sh).
threadNum="$THREAD_NUM_OVERRIDE"
export THREAD_NUM="$threadNum"

# Extract the name of the parent directory of $org1FASTA
org1DirName="$(basename "$(dirname "$org1FASTA")")"
org2DirName="$(basename "$(dirname "$org2FASTA")")"
org3DirName="$(basename "$(dirname "$org3FASTA")")"

org1FullName="${org1DirName}_1"
org2FullName="${org2DirName}_2"
org3FullName="${org3DirName}_3"

# make short names
org1ShortName=$(make_short_name "$org1FullName" "1")
org2ShortName=$(make_short_name "$org2FullName" "2")
org3ShortName=$(make_short_name "$org3FullName" "3")

if [ -n "${EVOSUBSTER_METADATA_INPUT:-}" ]; then
    mapfile -t display_names < <(python3 -c 'import json,sys; print("\n".join(e["short_name"] for e in json.load(open(sys.argv[1]))))' "$EVOSUBSTER_METADATA_INPUT")
    [ "${#display_names[@]}" -eq 3 ] || exit 1
    org1ShortName="${display_names[0]}"; org2ShortName="${display_names[1]}"; org3ShortName="${display_names[2]}"
fi

org1ID=$(extract_accession_from_path "$org1FASTA")
org2ID=$(extract_accession_from_path "$org2FASTA")
org3ID=$(extract_accession_from_path "$org3FASTA")

if [ -z "$org1ID" ] || [ -z "$org2ID" ] || [ -z "$org3ID" ]; then
    echo "Error: Could not extract accession IDs from FASTA paths." 1>&2
    exit 1
fi

# Cache keys must identify the exact selected assemblies. Downloader mode
# supplies them explicitly; standalone mode enables reuse only when every FASTA
# basename starts with a complete versioned GCA/GCF accession.
if [ "${#ACCESSION_IDS[@]}" -eq 3 ]; then
    cacheOrg1ID="${ACCESSION_IDS[0]}"
    cacheOrg2ID="${ACCESSION_IDS[1]}"
    cacheOrg3ID="${ACCESSION_IDS[2]}"
else
    cacheOrg1ID=$(extract_ncbi_accession_from_path "$org1FASTA")
    cacheOrg2ID=$(extract_ncbi_accession_from_path "$org2FASTA")
    cacheOrg3ID=$(extract_ncbi_accession_from_path "$org3FASTA")
fi

if [ -n "$TRAIN_CACHE_DIR" ] &&
   { [ -z "$cacheOrg1ID" ] || [ -z "$cacheOrg2ID" ] || [ -z "$cacheOrg3ID" ]; }; then
    echo "Train-cache disabled: exact versioned GCA/GCF accessions were not supplied and could not be extracted from every FASTA filename." >&2
fi

# Use config patterns to generate filenames
default_out_dir=$(get_dwl_config '.paths.out_dir')
if [ -z "$default_out_dir" ] || [ "$default_out_dir" = "null" ]; then
    echo "Error: .paths.out_dir is not set in dwl_config.yaml" >&2
    exit 1
fi
outDirBase="${OUT_DIR_OVERRIDE:-$default_out_dir}"
outDirBase=$(resolve_path "$outDirBase")
if [ ! -d "$outDirBase" ]; then
    echo "---making $outDirBase"
    if ! mkdir -p "$outDirBase"; then
        echo "Error: Unable to create output base directory $outDirBase" >&2
        exit 1
    fi
fi
MANIFEST_HELPER="$SCRIPT_DIR/lib/run_manifest.py"
manifest_args=("$outDirBase" "$DATE" "$threadNum"
    "$org1FASTA" "${gff_paths[0]}" "$cacheOrg1ID" "$org1ShortName"
    "$org2FASTA" "${gff_paths[1]}" "$cacheOrg2ID" "$org2ShortName"
    "$org3FASTA" "${gff_paths[2]}" "$cacheOrg3ID" "$org3ShortName")
# Serialize selection through manifest publication, including short-key collisions.
exec 8>"$outDirBase/.run-allocation.lock"
flock 8 || exit 1
run_path=$(python3 "$MANIFEST_HELPER" path "${manifest_args[@]}") || exit 1
mkdir -p "$run_path" || exit 1
exec 9>"$run_path/.run.lock"
flock -n 9 || { echo "Error: run is already active: $run_path" >&2; exit 1; }
python3 "$MANIFEST_HELPER" prepare "${manifest_args[@]}" --run-dir "$run_path" >/dev/null || exit 1
flock -u 8
exec 8>&-
outDirPath="$run_path"
# A resumed run retains its recorded display labels even if input folders moved.
mapfile -t run_names < <(python3 -c 'import json,sys; print("\n".join(e["short_name"] for e in json.load(open(sys.argv[1]))["organisms"]))' "$run_path/metadata/metadata_manifest.json")
org1ShortName="${run_names[0]}"; org2ShortName="${run_names[1]}"; org3ShortName="${run_names[2]}"


gcContent_org1=$(get_config '.patterns.gc_content' | sed "s/{org_short}/$org1ShortName/g" | sed "s/{date}/$DATE/g")
gcContent_org2=$(get_config '.patterns.gc_content' | sed "s/{org_short}/$org2ShortName/g" | sed "s/{date}/$DATE/g")
gcContent_org3=$(get_config '.patterns.gc_content' | sed "s/{org_short}/$org3ShortName/g" | sed "s/{date}/$DATE/g")

sbstRatio=$(get_config '.patterns.sbst_ratio' | sed "s/{date}/$DATE/g")

dbName="$org1ShortName""db_$DATE"
o2o12=$(get_config '.patterns.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g" | sed "s/{date}/$DATE/g")
o2o13=$(get_config '.patterns.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org3ShortName/g" | sed "s/{date}/$DATE/g")
train12=$(get_config '.patterns.train' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g" | sed "s/{date}/$DATE/g")
train13=$(get_config '.patterns.train' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org3ShortName/g" | sed "s/{date}/$DATE/g")

joinedFile=$(get_config '.patterns.joined' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g" | sed "s/{org3_short}/$org3ShortName/g" | sed "s/{date}/$DATE/g")
joinedFile_ncds=$(get_config '.patterns.joined_ncds' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g" | sed "s/{org3_short}/$org3ShortName/g" | sed "s/{date}/$DATE/g")

org2tsv="${org2ID}_${org2ShortName}_${DATE}.tsv"
org3tsv="${org3ID}_${org3ShortName}_${DATE}.tsv"
org2_dinuc_tsv="${org2ID}_${org2ShortName}_${DATE}_dinuc.tsv"
org3_dinuc_tsv="${org3ID}_${org3ShortName}_${DATE}_dinuc.tsv"
org2tsv_ncds="${org2ID}_${org2ShortName}_${DATE}_ncds.tsv"
org3tsv_ncds="${org3ID}_${org3ShortName}_${DATE}_ncds.tsv"
org2_dinuc_tsv_ncds="${org2ID}_${org2ShortName}_${DATE}_dinuc_ncds.tsv"
org3_dinuc_tsv_ncds="${org3ID}_${org3ShortName}_${DATE}_dinuc_ncds.tsv"

cd "$outDirPath" || exit 1
last_options_guard "$outDirPath/intermediateFiles" training || exit 1
REUSE_HELPER="$SCRIPT_DIR/lib/retry_reuse.py"
if [ -n "${EVOSUBSTER_REUSE_SOURCE:-}" ]; then
    [ -n "${EVOSUBSTER_REUSE_LOG:-}" ] || { echo "Error: EVOSUBSTER_REUSE_LOG is required." >&2; exit 1; }
    python3 "$REUSE_HELPER" import "$outDirPath" "$EVOSUBSTER_REUSE_SOURCE" "$EVOSUBSTER_REUSE_LOG" || exit 1
fi
echo "pwd: $(pwd)"

# Create subdirectory structure
imf="intermediateFiles"
figs_org2="figs/${org2ShortName}"
figs_org3="figs/${org3ShortName}"
stat_org2="statistics/${org2ShortName}"
stat_org3="statistics/${org3ShortName}"
stat_misc="statistics/misc"
mkdir -p \
    "$imf" \
    "$figs_org2/singlenuc/ratio" \
    "$figs_org2/singlenuc/log-ratio" \
    "$figs_org2/singlenuc/count" \
    "$figs_org2/dinuc" \
    "$figs_org3/singlenuc/ratio" \
    "$figs_org3/singlenuc/log-ratio" \
    "$figs_org3/singlenuc/count" \
    "$figs_org3/dinuc" \
    "$stat_org2/singlenuc" \
    "$stat_org2/dinuc" \
    "$stat_org3/singlenuc" \
    "$stat_org3/dinuc" \
    "$stat_misc"

# Prefix path variables with subdirectory locations
gcContent_org1="${stat_misc}/${gcContent_org1}"
gcContent_org2="${stat_misc}/${gcContent_org2}"
gcContent_org3="${stat_misc}/${gcContent_org3}"
sbstRatio="${stat_misc}/${sbstRatio}"
o2o12="${imf}/${o2o12}"
o2o13="${imf}/${o2o13}"
train12="${imf}/${train12}"
train13="${imf}/${train13}"
joinedFile="${imf}/${joinedFile}"
joinedFile_ncds="${imf}/${joinedFile_ncds}"
org2tsv="${stat_org2}/singlenuc/${org2tsv}"
org3tsv="${stat_org3}/singlenuc/${org3tsv}"
org2_dinuc_tsv="${stat_org2}/dinuc/${org2_dinuc_tsv}"
org3_dinuc_tsv="${stat_org3}/dinuc/${org3_dinuc_tsv}"
org2tsv_ncds="${stat_org2}/singlenuc/${org2tsv_ncds}"
org3tsv_ncds="${stat_org3}/singlenuc/${org3tsv_ncds}"
org2_dinuc_tsv_ncds="${stat_org2}/dinuc/${org2_dinuc_tsv_ncds}"
org3_dinuc_tsv_ncds="${stat_org3}/dinuc/${org3_dinuc_tsv_ncds}"

# GC content
echo "$(get_config '.messages.gc_content')"
if [ ! -e "$gcContent_org1" ]; then
echo "time bash $METRICS_DIR/gc_content.sh $org1FASTA >$gcContent_org1"
time bash "$METRICS_DIR/gc_content.sh" "$org1FASTA" >"$gcContent_org1" || exit 1
else
	echo "$gcContent_org1 already exists"
fi
if [ ! -e "$gcContent_org2" ]; then
echo "time bash $METRICS_DIR/gc_content.sh $org2FASTA >$gcContent_org2"
time bash "$METRICS_DIR/gc_content.sh" "$org2FASTA" >"$gcContent_org2" || exit 1
else
	echo "$gcContent_org2 already exists"
fi
if [ ! -e "$gcContent_org3" ]; then
echo "time bash $METRICS_DIR/gc_content.sh $org3FASTA >$gcContent_org3"
time bash "$METRICS_DIR/gc_content.sh" "$org3FASTA" >"$gcContent_org3" || exit 1
else
	echo "$gcContent_org3 already exists"
fi

# Coding region percentage (all three species, when annotations are available).
echo "$(get_config '.messages.coding_region_percentage')"
coding_pattern=$(get_config '.patterns.coding_region_percentage')
coding_short_names=("$org1ShortName" "$org2ShortName" "$org3ShortName")
for i in 0 1 2; do
    if [ "${gff_paths[$i]}" = "NO_GFF_FILE" ]; then
        echo "No GFF for ${coding_short_names[$i]}; skipping coding region percentage."
        continue
    fi
    coding_filename="${coding_pattern//\{org_short\}/${coding_short_names[$i]}}"
    coding_filename="${coding_filename//\{date\}/$DATE}"
    coding_output="$stat_misc/$coding_filename"
    if [ -s "$coding_output" ]; then
        echo "$coding_output already exists"
    else
        python3 "$METRICS_DIR/coding_region_percentage.py" \
            "${gff_fastas[$i]}" "${gff_paths[$i]}" --output "$coding_output" || exit 1
    fi
done

# Run last-train to check substitution percent identity between org1 and org2
echo "$(get_config '.options.checkIdt.enabled_message')"
if python3 "$REUSE_HELPER" pair-ready "$outDirPath" org1-org2; then
    echo "Using certified LAST 1648 training for org1-org2 (original provenance retained)."
else
    reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg1ID" "$cacheOrg2ID" "$train12" "$org1FASTA" "$org2FASTA" || exit 1
    time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org1FASTA" "$org2FASTA" "$org1ShortName" "$org2ShortName" "$imf" || exit 1
fi
# Run last-train to check substitution percent identity between org1 and org3
echo "$(get_config '.options.checkIdt.enabled_message')"
if python3 "$REUSE_HELPER" pair-ready "$outDirPath" org1-org3; then
    echo "Using certified LAST 1648 training for org1-org3 (original provenance retained)."
else
    reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg1ID" "$cacheOrg3ID" "$train13" "$org1FASTA" "$org3FASTA" || exit 1
    time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org1FASTA" "$org3FASTA" "$org1ShortName" "$org3ShortName" "$imf" || exit 1
fi
# Run last-train to check substitution percent identity between org2 and org3 (inner group)
echo "$(get_config '.options.checkIdt.enabled_message')"
reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg2ID" "$cacheOrg3ID" \
    "${imf}/${org2ShortName}2${org3ShortName}_${DATE}.train" "$org2FASTA" "$org3FASTA" || exit 1
time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org2FASTA" "$org3FASTA" "$org2ShortName" "$org3ShortName" "$imf" || exit 1

for train_file in "$train12" "$train13" "${imf}/${org2ShortName}2${org3ShortName}_${DATE}.train"; do
    is_valid_train_file "$train_file" || { echo "Error: incomplete training output: $train_file" >&2; exit 1; }
done
identity_specs=("gc:org1:$gcContent_org1" "gc:org2:$gcContent_org2" "gc:org3:$gcContent_org3")
for i in 0 1 2; do
    coding_filename="${coding_pattern//\{org_short\}/${coding_short_names[$i]}}"
    coding_filename="${coding_filename//\{date\}/$DATE}"
    if [ "${gff_paths[$i]}" != "NO_GFF_FILE" ]; then
        identity_specs+=("coding:org$((i+1)):$stat_misc/$coding_filename")
    fi
done
python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" identity_complete "${identity_specs[@]}" \
    "train:org1-org2:$train12" "train:org1-org3:$train13" \
    "train:org2-org3:${imf}/${org2ShortName}2${org3ShortName}_${DATE}.train" || exit 1

if [ "$IDT_ONLY" -eq 1 ]; then
    echo "---idt-only is enabled; exiting after last-train identity checks."
    exit 0
fi

# one2one for org1-org2
echo "$(get_config '.messages.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g")"
echo "bash $ALIGN_DIR/one2one.sh $DATE $org1FASTA $org2FASTA ${imf}/${dbName} $train12 $o2o12"
if python3 "$REUSE_HELPER" pair-ready "$outDirPath" org1-org2; then
    echo "Using certified LAST 1648 alignment for org1-org2 (original provenance retained)."
else
    last_options_guard "$o2o12" alignment || exit 1
    if python3 "$REUSE_HELPER" current-pair-ready "$outDirPath" org1-org2; then
        echo "Using checkpointed sorted alignment for org1-org2."
    else
        bash "$ALIGN_DIR/one2one.sh" "$DATE" "$org1FASTA" "$org2FASTA" "${imf}/${dbName}" "$train12" "$o2o12" || exit 1
    fi
fi
python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" identity_complete || exit 1

# one2one for org1-org3
echo "$(get_config '.messages.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org3ShortName/g")"
echo "bash $ALIGN_DIR/one2one.sh $DATE $org1FASTA $org3FASTA ${imf}/${dbName} $train13 $o2o13"
if python3 "$REUSE_HELPER" pair-ready "$outDirPath" org1-org3; then
    echo "Using certified LAST 1648 alignment for org1-org3 (original provenance retained)."
else
    last_options_guard "$o2o13" alignment || exit 1
    if python3 "$REUSE_HELPER" current-pair-ready "$outDirPath" org1-org3; then
        echo "Using checkpointed sorted alignment for org1-org3."
    else
        bash "$ALIGN_DIR/one2one.sh" "$DATE" "$org1FASTA" "$org3FASTA" "${imf}/${dbName}" "$train13" "$o2o13" || exit 1
    fi
fi
python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" identity_complete || exit 1

# maf-join the two .maf files
echo "$(get_config '.messages.maf_join')"
echo "bash $ALIGN_DIR/mafjoin.sh $o2o12 $o2o13 $joinedFile"
bash "$ALIGN_DIR/mafjoin.sh" "$o2o12" "$o2o13" "$joinedFile" || exit 1
python3 "$REUSE_HELPER" refresh-sorted "$outDirPath" || exit 1
python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" identity_complete || exit 1

# Calculate the substitution ratio without considering neighboring bases
echo "$(get_config '.messages.sbst_ratio')"
ratio_complete="${sbstRatio}.complete"
if [ ! -s "$sbstRatio" ] || [ ! -e "$ratio_complete" ]; then
echo "time python3 $METRICS_DIR/subRatio.py $joinedFile >$sbstRatio"
ratio_tmp="${sbstRatio}.tmp.$$"
if time python3 "$METRICS_DIR/subRatio.py" "$joinedFile" >"$ratio_tmp" \
        && [ -s "$ratio_tmp" ]; then
    mv "$ratio_tmp" "$sbstRatio"
    : > "$ratio_complete"
else
    rm -f "$ratio_tmp"
    exit 1
fi
else
    echo "$sbstRatio already exists"
fi

# If there is a gff file of org1, cut off the CDS regions
# and count the number of substitutions in non-coding regions
if [ "$org1GFF" != "NO_GFF_FILE" ]; then
	echo "There is a gff file of org1"
	echo "maf-cut (cut off the CDS regions)"
	if [ ! -e "$joinedFile_ncds" ]; then
		# Write to a tmp file first; rename on success so a failed run
		# (e.g. GFF/FASTA seq-name mismatch) leaves no partial ncds MAF
		# behind that the next run would blindly skip regenerating.
		tmp_ncds="${joinedFile_ncds}.tmp"
		if "$ALIGN_DIR/maf-cut-cds-uglier.py" \
				"$org1GFF" \
				"$joinedFile" >"$tmp_ncds"; then
			mv "$tmp_ncds" "$joinedFile_ncds"
		else
			rc=$?
			rm -f "$tmp_ncds"
			echo "Error: maf-cut-cds-uglier.py failed (exit $rc); aborting." >&2
			exit "$rc"
		fi
	else
		echo "$joinedFile_ncds already exists"
	fi

	# Generate all TSV files including ncds files
	bash "$SCRIPT_DIR/generate_tsv_files.sh" \
		"$joinedFile" \
		"$org2tsv" \
		"$org3tsv" \
		"$COUNT_DIR" \
		"$org2_dinuc_tsv" \
		"$org3_dinuc_tsv" \
		"$joinedFile_ncds" \
		"$org2tsv_ncds" \
		"$org3tsv_ncds" \
		"$org2_dinuc_tsv_ncds" \
		"$org3_dinuc_tsv_ncds" || exit 1
else
	echo "There is no gff file of org1"
	# Generate all TSV files (no ncds files)
	bash "$SCRIPT_DIR/generate_tsv_files.sh" \
		"$joinedFile" \
		"$org2tsv" \
		"$org3tsv" \
		"$COUNT_DIR" \
		"$org2_dinuc_tsv" \
		"$org3_dinuc_tsv" || exit 1
fi

count_specs=("gc:org1:$gcContent_org1" "gc:org2:$gcContent_org2" "gc:org3:$gcContent_org3" "sbst_ratio:all:$sbstRatio" "single:org2:$org2tsv" "single:org3:$org3tsv"
    "dinuc:org2:$org2_dinuc_tsv" "dinuc:org3:$org3_dinuc_tsv")
if [ "$org1GFF" != "NO_GFF_FILE" ]; then
    count_specs+=("single_ncds:org2:$org2tsv_ncds" "single_ncds:org3:$org3tsv_ncds"
        "dinuc_ncds:org2:$org2_dinuc_tsv_ncds" "dinuc_ncds:org3:$org3_dinuc_tsv_ncds")
fi
python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" counts_complete "${count_specs[@]}" || exit 1

# Generate all graphs
bash "$SCRIPT_DIR/generate_graphs.sh" \
    "$org2tsv" \
    "$org3tsv" \
    "$org2tsv_ncds" \
    "$org3tsv_ncds" \
    "$org2_dinuc_tsv" \
    "$org3_dinuc_tsv" \
    "$org2_dinuc_tsv_ncds" \
    "$org3_dinuc_tsv_ncds" \
    "$R_DIR" \
    "$figs_org2" \
    "$figs_org3" || exit 1

python3 "$MANIFEST_HELPER" checkpoint "$outDirPath" complete || exit 1
echo "[[ This is the end of the pipeline ]]"
