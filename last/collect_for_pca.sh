#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/sbst_config.yaml"
DWL_CONFIG_FILE="$SCRIPT_DIR/dwl_config.yaml"

usage() {
    echo "Usage: $0 <metadata_manifest.json> [--tsv-dir PATH] [--pca-root PATH]" >&2
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

manifestPath=$1
shift

tsvDirOverride=""
pcaRootOverride=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tsv-dir)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --tsv-dir requires a non-empty path argument." >&2
                exit 1
            fi
            tsvDirOverride="$2"
            shift 2
            ;;
        --tsv-dir=*)
            tsvDirOverride="${1#*=}"
            if [[ -z "$tsvDirOverride" ]]; then
                echo "Error: --tsv-dir requires a non-empty path argument." >&2
                exit 1
            fi
            shift
            ;;
        --pca-root)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --pca-root requires a non-empty path argument." >&2
                exit 1
            fi
            pcaRootOverride="$2"
            shift 2
            ;;
        --pca-root=*)
            pcaRootOverride="${1#*=}"
            if [[ -z "$pcaRootOverride" ]]; then
                echo "Error: --pca-root requires a non-empty path argument." >&2
                exit 1
            fi
            shift
            ;;
        *)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
    esac
done

if [ ! -f "$manifestPath" ]; then
    echo "Error: metadata manifest not found at $manifestPath" >&2
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: sbst_config.yaml not found at $CONFIG_FILE" >&2
    exit 1
fi

get_config() {
    yq eval "$1" "$CONFIG_FILE"
}

pcaRoot="$pcaRootOverride"
if [ -z "$pcaRoot" ]; then
    pcaRoot=$(get_config '.paths.pca_analysis')
fi

if [ -z "$pcaRoot" ] || [ "$pcaRoot" = "null" ]; then
    echo "Error: PCA root directory not provided and .paths.pca_analysis is not set in sbst_config.yaml" >&2
    exit 1
fi

if [[ "$pcaRoot" != /* ]]; then
    baseOutDir=$(get_config '.paths.out_dir')
    if [ -z "$baseOutDir" ] || [ "$baseOutDir" = "null" ]; then
        echo "Error: .paths.out_dir is not set in sbst_config.yaml" >&2
        exit 1
    fi
    pcaRoot="$baseOutDir/$pcaRoot"
fi

if [ ! -d "$pcaRoot" ]; then
    mkdir -p "$pcaRoot"
fi

pcaTsvDir="$pcaRoot/tsv"
pcaMetadataDir="$pcaRoot/metadata"
pcaResultsDir="$pcaRoot/results"

mkdir -p "$pcaTsvDir" "$pcaMetadataDir" "$pcaResultsDir"

manifestDate=$(jq -r '.date // empty' "$manifestPath")
manifestCombo=$(jq -r '.combo_dir // empty' "$manifestPath")
runDir=$(jq -r '.run_dir // empty' "$manifestPath")
metadataDir=$(jq -r '.metadata_dir // empty' "$manifestPath")

if [ -z "$manifestDate" ] || [ "$manifestDate" = "null" ]; then
    echo "Error: Unable to read 'date' from $manifestPath" >&2
    exit 1
fi

if [ -z "$metadataDir" ] || [ "$metadataDir" = "null" ]; then
    echo "Error: Unable to read 'metadata_dir' from $manifestPath" >&2
    exit 1
fi

if [ ! -d "$metadataDir" ]; then
    echo "Error: metadata_dir $metadataDir does not exist" >&2
    exit 1
fi

tsvSourceDir="$runDir"
if [ -n "$tsvDirOverride" ]; then
    tsvSourceDir="$tsvDirOverride"
fi

if [ -z "$tsvSourceDir" ]; then
    echo "Error: TSV directory could not be determined; specify --tsv-dir" >&2
    exit 1
fi

if [ ! -d "$tsvSourceDir" ]; then
    echo "Error: TSV directory $tsvSourceDir does not exist" >&2
    exit 1
fi

declare -a SHORT_NAMES=()
declare -a ACCESSIONS=()
declare -a METADATA_PATHS=()

while IFS=$'\t' read -r short accession metadata_json role; do
    if [ "$role" = "outgroup" ]; then
        continue
    fi
    SHORT_NAMES+=("$short")
    ACCESSIONS+=("$accession")
    METADATA_PATHS+=("$metadata_json")
done < <(jq -r '.organisms[] | [.short_name, .accession, (.metadata_json // ""), (.role // "")] | @tsv' "$manifestPath")

if [ ${#SHORT_NAMES[@]} -eq 0 ]; then
    echo "Error: No ingroup organisms found in $manifestPath" >&2
    exit 1
fi

link_tsv_files() {
    local accession=$1
    local short=$2
    local pattern="${accession}_${short}_${manifestDate}*.tsv"

    shopt -s nullglob
    local matches=("$tsvSourceDir"/$pattern)
    shopt -u nullglob

    if [ ${#matches[@]} -eq 0 ]; then
        echo "Warning: TSV files matching $pattern not found in $tsvSourceDir" >&2
        return
    fi

    for src in "${matches[@]}"; do
        local name
        name=$(basename "$src")
        ln -sf "$src" "$pcaTsvDir/$name"
    done
}

link_into_pca_metadata() {
    local src=$1
    if [ -z "$src" ] || [ ! -f "$src" ]; then
        return
    fi
    local name
    name=$(basename "$src")
    ln -sf "$src" "$pcaMetadataDir/$name"
}

ensure_taxonomy_json() {
    local short=$1
    local accession=$2
    local metadata_json=$3
    local taxonomy_dir=$4

    local taxonomy_file="$taxonomy_dir/taxonomy_${short}_${accession}.json"
    if [ -s "$taxonomy_file" ]; then
        echo "$taxonomy_file"
        return 0
    fi

    if [ -z "$metadata_json" ] || [ ! -f "$metadata_json" ]; then
        echo "Warning: Metadata JSON $metadata_json missing for ${short} (${accession}); cannot derive taxonomy" >&2
        return 1
    fi

    local tax_id
    tax_id=$(jq -r 'try .reports[0].organism.tax_id catch ""' "$metadata_json")
    if [ -z "$tax_id" ] || [ "$tax_id" = "null" ]; then
        tax_id=$(jq -r 'try .reports[0].biosample.description.organism.tax_id catch ""' "$metadata_json")
    fi

    if [ -z "$tax_id" ] || [ "$tax_id" = "null" ]; then
        echo "Warning: tax_id not found in $metadata_json for ${short} (${accession})" >&2
        return 1
    fi

    local tmp
    tmp=$(mktemp) || {
        echo "Warning: Unable to create temporary file for taxonomy download" >&2
        return 1
    }

    if datasets summary taxonomy taxon "$tax_id" >"$tmp"; then
        if mv "$tmp" "$taxonomy_file"; then
            echo "$taxonomy_file"
            return 0
        fi
    else
        echo "Warning: Failed to fetch taxonomy for ${short} (${accession}) with tax_id $tax_id" >&2
    fi

    rm -f "$tmp"
    return 1
}

for idx in "${!SHORT_NAMES[@]}"; do
    short=${SHORT_NAMES[$idx]}
    accession=${ACCESSIONS[$idx]}
    metadata_json=${METADATA_PATHS[$idx]}

    link_tsv_files "$accession" "$short"

    if [ -n "$metadata_json" ] && [ -f "$metadata_json" ]; then
        link_into_pca_metadata "$metadata_json"
    else
        echo "Warning: Metadata JSON missing for ${short} (${accession})" >&2
    fi

    taxonomy_json=$(ensure_taxonomy_json "$short" "$accession" "$metadata_json" "$metadataDir" || true)
    if [ -n "$taxonomy_json" ] && [ -f "$taxonomy_json" ]; then
        link_into_pca_metadata "$taxonomy_json"
    fi
done

if [ -n "$manifestCombo" ] && [ "$manifestCombo" != "null" ]; then
    ln -sf "$manifestPath" "$pcaMetadataDir/manifest_${manifestCombo}_${manifestDate}.json"
else
    ln -sf "$manifestPath" "$pcaMetadataDir/manifest_${manifestDate}.json"
fi

# Optionally run PCA automation if script exists and environment variable is set
PCA_AUTO_SCRIPT="$SCRIPT_DIR/../analysis/pca_auto.py"
if [ "${RUN_PCA_AUTO:-0}" -eq 1 ] && [ -f "$PCA_AUTO_SCRIPT" ]; then
    echo "--- running pca_auto.py"
    python "$PCA_AUTO_SCRIPT" \
        --tsv_dir "$pcaTsvDir" \
        --metadata_dir "$pcaMetadataDir" \
        --output_dir "$pcaResultsDir"
fi
