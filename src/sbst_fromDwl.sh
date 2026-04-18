#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PATH="$SCRIPT_DIR:$PATH"

config_file="$ROOT_DIR/config/dwl_config.yaml"
# Load YAML configuration using yq
if [ ! -f "$config_file" ]; then
    echo "Configuration file not found!" 1>&2
    exit 1
fi
# Function to get config values using yq
get_config() {
    yq eval "$1" "$config_file"
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


OUT_DIR_OVERRIDE=""
BASE_GENOMES_OVERRIDE=""
IDT_ONLY=0
THREAD_NUM_OVERRIDE=8
FORCE_DOWNLOAD=0
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --idt-only)
            IDT_ONLY=1
            shift
            continue
            ;;
        --force)
            FORCE_DOWNLOAD=1
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
        --genomes-dir)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --genomes-dir requires a non-empty path argument." >&2
                exit 1
            fi
            BASE_GENOMES_OVERRIDE="$2"
            shift 2
            continue
            ;;
        --genomes-dir=*)
            BASE_GENOMES_OVERRIDE="${1#*=}"
            if [[ -z "$BASE_GENOMES_OVERRIDE" ]]; then
                echo "Error: --genomes-dir requires a non-empty path argument." >&2
                exit 1
            fi
            shift
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

# Parse positional arguments
DATE="$1"
org1ID="$2"
org2ID="$3"
org3ID="$4"


# Check minimally required arguments (DATE and 3 accessions)
if [ -z "$DATE" ] || [ -z "$org1ID" ] || [ -z "$org2ID" ] || [ -z "$org3ID" ]; then
    echo "$(get_config '.errors.arg_count' | sed "s/{arg_num}/$(get_config '.settings.required_args')/g")" >&2
   echo "$(get_config '.errors.usage')" >&2
   exit 1
fi

default_base_genomes=$(get_config '.paths.base_genomes')
if [ -z "$BASE_GENOMES_OVERRIDE" ] && ([ -z "$default_base_genomes" ] || [ "$default_base_genomes" = "null" ]); then
    echo "Error: .paths.base_genomes is not set in dwl_config.yaml" >&2
    exit 1
fi
base_genomes="${BASE_GENOMES_OVERRIDE:-$default_base_genomes}"
base_genomes=$(resolve_path "$base_genomes")

if [ ! -d "$base_genomes" ]; then
    mkdir -p "$base_genomes" || {
        echo "Error: Unable to create base genomes directory at $base_genomes" >&2
        exit 1
    }
fi

default_out_dir=$(get_config '.paths.out_dir')
if [ -z "$default_out_dir" ] || [ "$default_out_dir" = "null" ]; then
    echo "Error: .paths.out_dir is not set in dwl_config.yaml" >&2
    exit 1
fi

out_dir_base="${OUT_DIR_OVERRIDE:-$default_out_dir}"
out_dir_base=$(resolve_path "$out_dir_base")
if [ ! -d "$out_dir_base" ]; then
    if ! mkdir -p "$out_dir_base"; then
        echo "Error: Unable to create output base directory at $out_dir_base" >&2
        exit 1
    fi
fi

# Download genome + taxonomy for each accession via dwl_organism.sh
_dwl_flags=()
[ "$FORCE_DOWNLOAD" -eq 1 ] && _dwl_flags+=("--force")
org1Result=$(bash "$SCRIPT_DIR/dwl_organism.sh" "$org1ID" --out-dir "$base_genomes" "${_dwl_flags[@]}") || exit 1
org2Result=$(bash "$SCRIPT_DIR/dwl_organism.sh" "$org2ID" --out-dir "$base_genomes" "${_dwl_flags[@]}") || exit 1
org3Result=$(bash "$SCRIPT_DIR/dwl_organism.sh" "$org3ID" --out-dir "$base_genomes" "${_dwl_flags[@]}") || exit 1

IFS='|' read -r org1FullName org1FASTA org1SummaryJson org1RawName org1NcbiFullName org1TaxJson <<< "$org1Result"
IFS='|' read -r org2FullName org2FASTA org2SummaryJson org2RawName org2NcbiFullName org2TaxJson <<< "$org2Result"
IFS='|' read -r org3FullName org3FASTA org3SummaryJson org3RawName org3NcbiFullName org3TaxJson <<< "$org3Result"

echo "Derived org1FullName: $org1FullName"
echo "Derived org2FullName: $org2FullName"
echo "Derived org3FullName: $org3FullName"

# Create arrays
ids=("$org1ID" "$org2ID" "$org3ID")
names=("$org1FullName" "$org2FullName" "$org3FullName")
raw_names=("$org1RawName" "$org2RawName" "$org3RawName")
ncbi_full_names=("$org1NcbiFullName" "$org2NcbiFullName" "$org3NcbiFullName")
summary_jsons=("$org1SummaryJson" "$org2SummaryJson" "$org3SummaryJson")
tax_jsons=("$org1TaxJson" "$org2TaxJson" "$org3TaxJson")
fasta_paths=("$org1FASTA" "$org2FASTA" "$org3FASTA")

for f in "$org1FASTA" "$org2FASTA" "$org3FASTA"; do
    if [ -z "$f" ]; then
        echo "Error: Could not determine FASTA file path(s). Please ensure .fna files exist in the respective genome directories." >&2
        exit 1
    fi
done

org1FullNameWithSuffix="${org1FullName}_1"
org2FullNameWithSuffix="${org2FullName}_2"
org3FullNameWithSuffix="${org3FullName}_3"

org1ShortName=$(make_short_name "$org1FullNameWithSuffix" "1")
org2ShortName=$(make_short_name "$org2FullNameWithSuffix" "2")
org3ShortName=$(make_short_name "$org3FullNameWithSuffix" "3")

short_names=("$org1ShortName" "$org2ShortName" "$org3ShortName")
run_combo_dir="${org1ShortName}_${org2ShortName}_${org3ShortName}"
run_dir_root="$out_dir_base/$run_combo_dir"
run_date_dir_rel="$run_dir_root/$DATE"
metadata_dir_rel="$run_date_dir_rel/metadata"

if ! mkdir -p "$metadata_dir_rel"; then
    echo "Error: Unable to create metadata directory at $metadata_dir_rel" >&2
    exit 1
fi

run_date_dir="$(cd "$run_dir_root/$DATE" && pwd)"
metadata_dir="$run_date_dir/metadata"

# Check if GFF file exists in org1 (auto-detect)
gffFilePath="$base_genomes/$org1FullName/genomic.gff"

if [ -e "$gffFilePath" ]; then
    org1GFF="$gffFilePath" # set $org1GFF as the path to the gff file
    echo "GFF file found: $org1GFF"
else
    org1GFF="NO_GFF_FILE" # set a special flag to $org1GFF
    echo "No GFF file found for $org1FullName. Please download the GFF file manually."
fi

echo "--- Preparing metadata artifacts in $metadata_dir"
declare -a gff_paths=("$org1GFF" "" "")
declare -a metadata_json_copies=("" "" "")

for i in {0..2}; do
    src_json=${summary_jsons[$i]}
    short_name=${short_names[$i]}
    accession=${ids[$i]}
    if [ -n "$src_json" ] && [ -f "$src_json" ]; then
        dest_json="$metadata_dir/${short_name}_${accession}.json"
        if cp "$src_json" "$dest_json"; then
            metadata_json_copies[$i]="$dest_json"
        else
            echo "Warning: Failed to copy metadata JSON for ${short_name} ($accession) to $dest_json" >&2
        fi
    else
        echo "Warning: Source metadata JSON missing for ${short_name} ($accession) at ${src_json:-<empty>}" >&2
    fi
done

manifest_path="$metadata_dir/metadata_manifest.json"
python - "$manifest_path" "$DATE" "$run_combo_dir" "$run_date_dir" "$metadata_dir" \
"${ids[0]}" "${names[0]}" "${short_names[0]}" "${fasta_paths[0]}" "${gff_paths[0]}" "${raw_names[0]}" "${ncbi_full_names[0]}" "${metadata_json_copies[0]}" \
"${ids[1]}" "${names[1]}" "${short_names[1]}" "${fasta_paths[1]}" "${gff_paths[1]}" "${raw_names[1]}" "${ncbi_full_names[1]}" "${metadata_json_copies[1]}" \
"${ids[2]}" "${names[2]}" "${short_names[2]}" "${fasta_paths[2]}" "${gff_paths[2]}" "${raw_names[2]}" "${ncbi_full_names[2]}" "${metadata_json_copies[2]}" <<'PY'
import json
import sys
from pathlib import Path

def normalize(value, empty_as_none=False):
    if value in ("", "null", "None"):
        return None
    if empty_as_none and value in ("NO_GFF_FILE",):
        return None
    return value

manifest_path = Path(sys.argv[1])
run_date = sys.argv[2]
combo_dir = sys.argv[3]
run_dir = sys.argv[4]
metadata_dir = sys.argv[5]

args = sys.argv[6:]
slots = [
    ("org1", "outgroup"),
    ("org2", "ingroup"),
    ("org3", "ingroup"),
]
chunk_size = 8
organisms = []
for idx, (slot, role) in enumerate(slots):
    offset = idx * chunk_size
    chunk = args[offset : offset + chunk_size]
    accession, directory_name, short_name, fasta_path, gff_path, raw_name, ncbi_full, metadata_json = chunk
    organisms.append(
        {
            "slot": slot,
            "role": role,
            "accession": normalize(accession),
            "directory_name": normalize(directory_name),
            "short_name": normalize(short_name),
            "fasta_path": normalize(fasta_path),
            "gff_path": normalize(gff_path, empty_as_none=True),
            "raw_organism_name": normalize(raw_name),
            "ncbi_full_name": normalize(ncbi_full),
            "metadata_json": normalize(metadata_json),
        }
    )

manifest = {
    "date": run_date,
    "combo_dir": combo_dir,
    "run_dir": run_dir,
    "metadata_dir": metadata_dir,
    "organisms": organisms,
}
manifest_path.parent.mkdir(parents=True, exist_ok=True)
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"Wrote metadata manifest to {manifest_path}")
PY

# shellcheck disable=SC2181
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate metadata manifest" >&2
    exit 1
fi

# Copy taxonomy JSONs (downloaded by dwl_organism.sh) to metadata dir with short_name prefix
for i in {0..2}; do
    tax_src=${tax_jsons[$i]}
    short_name=${short_names[$i]}
    accession=${ids[$i]}
    tax_dest="$metadata_dir/taxonomy_${short_name}_${accession}.json"
    if [ -f "$tax_src" ] && [ ! -f "$tax_dest" ]; then
        cp "$tax_src" "$tax_dest"
    fi
done

# Run downstream pipeline (no checkInnerGroupIdt argument anymore)
trisbst_args=()
if [ -n "$OUT_DIR_OVERRIDE" ]; then
    trisbst_args+=("--out-dir" "$OUT_DIR_OVERRIDE")
fi
trisbst_args+=("--thread" "$THREAD_NUM_OVERRIDE")
if [ "$IDT_ONLY" -eq 1 ]; then
    trisbst_args+=("--idt-only")
fi
trisbst_args+=("$DATE" "$org1FASTA" "$org2FASTA" "$org3FASTA" "$org1GFF")
bash "$SCRIPT_DIR/sbst.sh" "${trisbst_args[@]}"
