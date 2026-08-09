#!/bin/bash

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
  --train-cache-dir PATH      Directory containing cached last-train files
  --accession-ids A1 A2 A3   Exact NCBI accessions for cache keys; normally supplied by sbst_fromDwl.sh
  --thread N                  Number of threads for LAST alignment (default: 8)
  --idt-only                  Stop after checking sequence percent identity among three genomes; skip downstream analysis
  -h, --help                  Show this help message and exit
EOF
        exit 0
    fi
done

lastal --version

# Resolve script locations and configuration relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
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
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
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

# Propagate thread count to child scripts (last_train.sh / one2one.sh).
threadNum="$THREAD_NUM_OVERRIDE"
export THREAD_NUM="$threadNum"

# Extract the name of the parent directory of $org1FASTA
org1DirName="$(basename $(dirname $org1FASTA))"
org2DirName="$(basename $(dirname $org2FASTA))"
org3DirName="$(basename $(dirname $org3FASTA))"

org1FullName="${org1DirName}_1"
org2FullName="${org2DirName}_2"
org3FullName="${org3DirName}_3"

# make short names
org1ShortName=$(make_short_name "$org1FullName" "1")
org2ShortName=$(make_short_name "$org2FullName" "2")
org3ShortName=$(make_short_name "$org3FullName" "3")

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
outDirPath="${outDirBase}/${org1ShortName}_${org2ShortName}_${org3ShortName}"

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

if [ ! -d "$outDirPath" ]; then
	echo "---making $outDirPath"
	mkdir -p "$outDirPath"
fi
cd "$outDirPath"
if [ ! -d "$DATE" ]; then
	echo "---making $DATE"
	mkdir "$DATE"
fi
cd "$DATE"
outDirPath=$(pwd)
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
time bash "$METRICS_DIR/gc_content.sh" "$org1FASTA" >"$gcContent_org1"
else
	echo "$gcContent_org1 already exists"
fi
if [ ! -e "$gcContent_org2" ]; then
echo "time bash $METRICS_DIR/gc_content.sh $org2FASTA >$gcContent_org2"
time bash "$METRICS_DIR/gc_content.sh" "$org2FASTA" >"$gcContent_org2"
else
	echo "$gcContent_org2 already exists"
fi
if [ ! -e "$gcContent_org3" ]; then
echo "time bash $METRICS_DIR/gc_content.sh $org3FASTA >$gcContent_org3"
time bash "$METRICS_DIR/gc_content.sh" "$org3FASTA" >"$gcContent_org3"
else
	echo "$gcContent_org3 already exists"
fi

# Run last-train to check substitution percent identity between org1 and org2
echo "$(get_config '.options.checkIdt.enabled_message')"
reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg1ID" "$cacheOrg2ID" "$train12" || exit 1
time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org1FASTA" "$org2FASTA" "$org1ShortName" "$org2ShortName" "$imf"
# Run last-train to check substitution percent identity between org1 and org3
echo "$(get_config '.options.checkIdt.enabled_message')"
reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg1ID" "$cacheOrg3ID" "$train13" || exit 1
time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org1FASTA" "$org3FASTA" "$org1ShortName" "$org3ShortName" "$imf"
# Run last-train to check substitution percent identity between org2 and org3 (inner group)
echo "$(get_config '.options.checkIdt.enabled_message')"
reuse_cached_train "$TRAIN_CACHE_DIR" "$DATE" "$cacheOrg2ID" "$cacheOrg3ID" \
    "${imf}/${org2ShortName}2${org3ShortName}_${DATE}.train" || exit 1
time bash "$ALIGN_DIR/last_train.sh" "$DATE" "$org2FASTA" "$org3FASTA" "$org2ShortName" "$org3ShortName" "$imf"

if [ "$IDT_ONLY" -eq 1 ]; then
    echo "---idt-only is enabled; exiting after last-train identity checks."
    exit 0
fi

# one2one for org1-org2
echo "$(get_config '.messages.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org2ShortName/g")"
echo "bash $ALIGN_DIR/one2one.sh $DATE $org1FASTA $org2FASTA ${imf}/${dbName} $train12 $o2o12"
bash "$ALIGN_DIR/one2one.sh" "$DATE" "$org1FASTA" "$org2FASTA" "${imf}/${dbName}" "$train12" "$o2o12" || exit 1

# one2one for org1-org3
echo "$(get_config '.messages.one2one' | sed "s/{org1_short}/$org1ShortName/g" | sed "s/{org2_short}/$org3ShortName/g")"
echo "bash $ALIGN_DIR/one2one.sh $DATE $org1FASTA $org3FASTA ${imf}/${dbName} $train13 $o2o13"
bash "$ALIGN_DIR/one2one.sh" "$DATE" "$org1FASTA" "$org3FASTA" "${imf}/${dbName}" "$train13" "$o2o13" || exit 1

# maf-join the two .maf files
echo "$(get_config '.messages.maf_join')"
echo "bash $ALIGN_DIR/mafjoin.sh $o2o12 $o2o13 $joinedFile"
bash "$ALIGN_DIR/mafjoin.sh" "$o2o12" "$o2o13" "$joinedFile" || exit 1

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

echo "[[ This is the end of the pipeline ]]"
