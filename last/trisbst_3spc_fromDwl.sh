#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LAST_DIR="$SCRIPT_DIR"

LAST_DIR="${LAST_DIR_OVERRIDE:-$LAST_DIR}"

PATH="$LAST_DIR:$PATH"

config_file="$SCRIPT_DIR/dwl_config.yaml"
# Load YAML configuration using yq
if [ ! -f "$config_file" ]; then
    echo "Configuration file not found!" 1>&2
    exit 1
fi
# Function to get config values using yq
get_config() {
    yq eval "$1" "$config_file"
}

sanitize_for_path() {
    local name="$1"
    name=${name// /_}
    name=${name//[^A-Za-z0-9._-]/_}
    echo "$name"
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

# Function to derive organism full name from NCBI Datasets summary JSON
# - Writes summary JSON to "$base_genomes/<orgFullName>/<accession>.json" (final_dir_name honors overrides)
# - Returns: "<final_dir_name>|<summary_json_path>|<raw_organism_name>|<calculated_full_name>"
get_org_full_name_from_id() {
    local accession="$1"
    local override_name="$2"
    local tmp_json
    tmp_json=$(mktemp) || {
        echo "Error: Unable to create temporary file for $accession summary." >&2
        return 1
    }

    # Require jq
    if ! command -v jq >/dev/null 2>&1; then
        echo "Error: 'jq' is required but not found in PATH." >&2
        rm -f "$tmp_json"
        return 1
    fi

    # Fetch summary JSON
    if ! datasets summary genome accession "$accession" > "$tmp_json"; then
        echo "Error: Failed to run 'datasets summary' for $accession" >&2
        rm -f "$tmp_json"
        return 1
    fi

    # Parse organism name
    local raw_base_name
    raw_base_name=$(jq -r 'try .reports[0].organism.organism_name catch ""' "$tmp_json")
    if [ -z "$raw_base_name" ] || [ "$raw_base_name" = "null" ]; then
        echo "Warning: '.reports[0].organism.organism_name' not found in temporary summary; using accession $accession" >&2
        raw_base_name="$accession"
    fi
    local base_name
    base_name=$(sanitize_for_path "$raw_base_name")

    # Parse infraspecific names (optional)
    local infra
    infra=$(jq -r 'try ([.reports[0].organism.infraspecific_names[]] | map(tostring) | join("_")) catch ""' "$tmp_json")
    infra=$(sanitize_for_path "$infra")

    local calculated_full_name
    if [ -n "$infra" ] && [ "$infra" != "null" ]; then
        calculated_full_name="${base_name}_${infra}"
    else
        calculated_full_name="$base_name"
    fi
    calculated_full_name=$(sanitize_for_path "$calculated_full_name")

    local final_dir_name
    if [ -n "$override_name" ]; then
        final_dir_name=$(sanitize_for_path "$override_name")
    else
        final_dir_name="$calculated_full_name"
    fi

    # Determine final JSON location
    local destination_json=""
    if [ -n "$base_genomes" ] && [ "$base_genomes" != "null" ]; then
        local dest_dir="$base_genomes/$final_dir_name"
        if ! mkdir -p "$dest_dir"; then
            echo "Error: Unable to create directory $dest_dir for storing summary JSON." >&2
            rm -f "$tmp_json"
            return 1
        fi
        destination_json="$dest_dir/${accession}.json"
    else
        destination_json="${accession}.json"
    fi

    if ! mv "$tmp_json" "$destination_json"; then
        echo "Warning: Failed to move summary JSON to $destination_json; attempting to copy instead." >&2
        if ! cp "$tmp_json" "$destination_json"; then
            echo "Error: Unable to place summary JSON for $accession." >&2
            rm -f "$tmp_json"
            return 1
        fi
        rm -f "$tmp_json"
    fi

    echo "$final_dir_name|$destination_json|$raw_base_name|$calculated_full_name"
}

OUT_DIR_OVERRIDE=""
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
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

base_genomes=$(get_config '.paths.base_genomes')
if [ -z "$base_genomes" ] || [ "$base_genomes" = "null" ]; then
    echo "Error: .paths.base_genomes is not set in dwl_config.yaml" >&2
    exit 1
fi

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
if [ ! -d "$out_dir_base" ]; then
    if ! mkdir -p "$out_dir_base"; then
        echo "Error: Unable to create output base directory at $out_dir_base" >&2
        exit 1
    fi
fi

# Auto-generate org full names from NCBI Datasets summary (always capture metadata)
org1Metadata=$(get_org_full_name_from_id "$org1ID" "") || exit 1
org2Metadata=$(get_org_full_name_from_id "$org2ID" "") || exit 1
org3Metadata=$(get_org_full_name_from_id "$org3ID" "") || exit 1

IFS='|' read -r org1FullName org1SummaryJson org1RawName org1NcbiFullName <<< "$org1Metadata"
IFS='|' read -r org2FullName org2SummaryJson org2RawName org2NcbiFullName <<< "$org2Metadata"
IFS='|' read -r org3FullName org3SummaryJson org3RawName org3NcbiFullName <<< "$org3Metadata"

echo "Derived org1FullName: $org1FullName"
echo "Derived org2FullName: $org2FullName"
echo "Derived org3FullName: $org3FullName"

cd "$base_genomes" || {
    echo "Error: Cannot change directory to $base_genomes" >&2
    exit 1
}
for orgFullName in $org1FullName $org2FullName $org3FullName; do
    if [ ! -d "$orgFullName" ]; then
        mkdir "$orgFullName"
    fi
done

includes=$(get_config '.download.includes'|tr -d ' '|tr '\n' ','|sed 's/,$//')
echo "includes: $includes"

# Create arrays
ids=("$org1ID" "$org2ID" "$org3ID")
names=("$org1FullName" "$org2FullName" "$org3FullName")
raw_names=("$org1RawName" "$org2RawName" "$org3RawName")
ncbi_full_names=("$org1NcbiFullName" "$org2NcbiFullName" "$org3NcbiFullName")
summary_jsons=("$org1SummaryJson" "$org2SummaryJson" "$org3SummaryJson")

# Iterate over both arrays using an index
for i in {0..2}; do
    orgID=${ids[$i]}
    orgFullName=${names[$i]}
    if [ ! -e "$base_genomes/$orgFullName/ncbi_dataset.zip" ]; then
        echo "$(get_config '.messages.download' | sed "s/{org_full}/$orgFullName/g")"
        cd "$base_genomes/$orgFullName"
        
        # Run download and check if it succeeded
        if ! datasets download genome accession "$orgID" --include "$includes"; then
            echo "Error: Failed to download genome for $orgFullName (ID: $orgID)" >&2
            echo "Exiting process..." >&2
            exit 1
        fi
        
        # Verify the downloaded file exists and is not empty
        if [ ! -s "ncbi_dataset.zip" ]; then
            echo "Error: Download completed but ncbi_dataset.zip is empty or missing for $orgFullName" >&2
            echo "Exiting process..." >&2
            exit 1
        fi
        
        echo "Successfully downloaded genome for $orgFullName"
    else
        echo "$(get_config '.messages.already_downloaded' | sed "s/{org_full}/$orgFullName/g")"
    fi
done

echo "All downloads completed successfully"

# move files and delete unnecessary directories
echo "$(get_config '.messages.move_files')"
function processGenomeData() {
    local orgFullName=$1
    local orgID=$2

    cd "$(get_config '.paths.base_genomes')/$orgFullName" || {
        echo "Error: Cannot change directory to $orgFullName" >&2
        return 1
    }

    if [ -z "$(ls *.fna 2>/dev/null)" ]; then
        echo "Processing genome data for $orgFullName..."
        
        # Unzip with error checking
        if ! unzip ncbi_dataset.zip; then
            echo "Error: Failed to unzip data for $orgFullName" >&2
            return 1
        fi

        cd ncbi_dataset/data || {
            echo "Error: Cannot access data directory for $orgFullName" >&2
            return 1
        }

        # Move files with error checking
        mv $(ls -p | grep -v /) "$(get_config '.paths.base_genomes')/$orgFullName" || {
            echo "Error: Failed to move files for $orgFullName" >&2
            return 1
        }

        cd "$orgID" || {
            echo "Error: Cannot access $orgID directory" >&2
            return 1
        }

        mv * "$(get_config '.paths.base_genomes')/$orgFullName" || {
            echo "Error: Failed to move $orgID files" >&2
            return 1
        }

        cd "$(get_config '.paths.base_genomes')/$orgFullName" || return 1
        rm -r ncbi_dataset || {
            echo "Warning: Could not remove ncbi_dataset directory" >&2
        }

        echo "Successfully processed genome data for $orgFullName"
    else
        echo "Genome files already exist for $orgFullName"
    fi
}

# Process each genome sequentially
for i in {0..2}; do
    orgFullName=${names[$i]}
    orgID=${ids[$i]}
    
    echo "Processing genome $((i+1)) of 3: $orgFullName"
    if ! processGenomeData "$orgFullName" "$orgID"; then
        echo "Error: Failed to process genome data for $orgFullName" >&2
        echo "Exiting process..." >&2
        exit 1
    fi
done

echo "All genome data processed successfully"

fasta_pat1=$(get_config '.patterns.fasta' | sed "s/{org_id}/$org1ID/g")
fasta_pat2=$(get_config '.patterns.fasta' | sed "s/{org_id}/$org2ID/g")
fasta_pat3=$(get_config '.patterns.fasta' | sed "s/{org_id}/$org3ID/g")

org1FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org1FullName/$fasta_pat1" | head -n1)
org2FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org2FullName/$fasta_pat2" | head -n1)
org3FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org3FullName/$fasta_pat3" | head -n1)

# Fallback to first *.fna if accession-specific pattern not found
if [ -z "$org1FASTA" ]; then org1FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org1FullName/*.fna" | head -n1); fi
if [ -z "$org2FASTA" ]; then org2FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org2FullName/*.fna" | head -n1); fi
if [ -z "$org3FASTA" ]; then org3FASTA=$(compgen -G "$(get_config '.paths.base_genomes')/$org3FullName/*.fna" | head -n1); fi

echo "org1FASTA: ${org1FASTA:-NOT_FOUND}"
echo "org2FASTA: ${org2FASTA:-NOT_FOUND}"
echo "org3FASTA: ${org3FASTA:-NOT_FOUND}"

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
gffFilePath="$(get_config '.paths.base_genomes')/$org1FullName/genomic.gff"

if [ -e "$gffFilePath" ]; then
    org1GFF="$gffFilePath" # set $org1GFF as the path to the gff file
    echo "GFF file found: $org1GFF"
else
    org1GFF="NO_GFF_FILE" # set a special flag to $org1GFF
    echo "No GFF file found for $org1FullName. Please download the GFF file manually."
fi

echo "--- Preparing metadata artifacts in $metadata_dir"
declare -a fasta_paths=("$org1FASTA" "$org2FASTA" "$org3FASTA")
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

# Run downstream pipeline (no checkInnerGroupIdt argument anymore)
trisbst_args=()
if [ -n "$OUT_DIR_OVERRIDE" ]; then
    trisbst_args+=("--out-dir" "$OUT_DIR_OVERRIDE")
fi
trisbst_args+=("$DATE" "$org1FASTA" "$org2FASTA" "$org3FASTA" "$org1GFF")
bash "$LAST_DIR/trisbst_3spc.sh" "${trisbst_args[@]}"
