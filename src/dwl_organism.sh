#!/bin/bash
# Download genome and taxonomy for a single NCBI accession.
#
# Usage: dwl_organism.sh <ACCESSION> [OPTIONS]
#
# Options:
#   --out-dir DIR    Base directory for genome storage (default: current directory)
#   --name NAME      Override organism directory name (default: derived from NCBI)
#   --no-genome      Skip genome download
#   --no-taxonomy    Skip taxonomy download
#   --tax-out FILE   Taxonomy JSON output path (default: <out-dir>/<dir_name>/taxonomy_<accession>.json)
#   --artifact-manifest FILE
#                    Write newly created artifacts as an atomic TSV manifest
#
# Stdout on success (pipe-separated):
#   dir_name|fasta_path|summary_json_path|raw_organism_name|ncbi_full_name|taxonomy_json_path

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- helpers ---

sanitize_for_path() {
    local name="$1"
    name=${name// /_}
    name=${name//[^A-Za-z0-9._-]/_}
    echo "$name"
}

curl_with_retry() {
    local attempt=1
    local max_attempts=4
    local status
    local retry_option=""

    while true; do
        if [ -n "$retry_option" ]; then
            curl "$retry_option" "$@" && return 0
            status=$?
        elif curl "$@"; then
            return 0
        else
            status=$?
        fi

        case "$status" in
            5|6|7|18|22|28|35|52|55|56|92) ;;
            *) return "$status" ;;
        esac
        [ "$attempt" -ge "$max_attempts" ] && return "$status"

        echo "Warning: NCBI request failed (curl exit $status); retrying with HTTP/1.1 ($((attempt + 1))/$max_attempts)..." >&2
        sleep 2
        retry_option="--http1.1"
        attempt=$((attempt + 1))
    done
}

snapshot_relative_paths() {
    local root="$1"
    local output="$2"

    : > "$output" || return 1
    if [ -d "$root" ]; then
        (
            cd "$root" || exit 1
            find . -mindepth 1 -print | LC_ALL=C sort
        ) > "$output"
    fi
}

write_artifact_manifest() {
    local manifest_path="$1"
    local org_dir="$2"
    local dir_name="$3"
    local org_dir_preexisting="$4"
    local before_snapshot="$5"
    local manifest_dir
    local path
    local relative_path
    local artifact_type

    ARTIFACT_AFTER_SNAPSHOT=$(mktemp) || {
        echo "Error: Unable to create artifact snapshot." >&2
        return 1
    }
    snapshot_relative_paths "$org_dir" "$ARTIFACT_AFTER_SNAPSHOT" || {
        echo "Error: Unable to snapshot artifacts in $org_dir" >&2
        return 1
    }

    manifest_dir=$(dirname "$manifest_path")
    mkdir -p "$manifest_dir" || {
        echo "Error: Unable to create artifact manifest directory $manifest_dir" >&2
        return 1
    }
    ARTIFACT_TMP_MANIFEST=$(mktemp "${manifest_path}.tmp.XXXXXX") || {
        echo "Error: Unable to create temporary artifact manifest for $manifest_path" >&2
        return 1
    }

    printf 'accession\tartifact_type\trelative_path\n' > "$ARTIFACT_TMP_MANIFEST" || {
        echo "Error: Unable to write artifact manifest $manifest_path" >&2
        return 1
    }
    if [ "$org_dir_preexisting" -eq 0 ]; then
        printf '%s\tdirectory_tree\t%s\n' "$ACCESSION" "$dir_name" \
            >> "$ARTIFACT_TMP_MANIFEST" || {
            echo "Error: Unable to write artifact manifest $manifest_path" >&2
            return 1
        }
    else
        while IFS= read -r path; do
            [ -n "$path" ] || continue
            if grep -Fqx -- "$path" "$before_snapshot"; then
                continue
            fi

            relative_path=${path#./}
            if [ -L "$org_dir/$relative_path" ]; then
                artifact_type="symlink"
            elif [ -d "$org_dir/$relative_path" ]; then
                artifact_type="directory"
            elif [ -f "$org_dir/$relative_path" ]; then
                artifact_type="file"
            else
                echo "Error: Unsupported artifact type at $org_dir/$relative_path" >&2
                return 1
            fi
            printf '%s\t%s\t%s/%s\n' \
                "$ACCESSION" "$artifact_type" "$dir_name" "$relative_path" \
                >> "$ARTIFACT_TMP_MANIFEST" || {
                echo "Error: Unable to write artifact manifest $manifest_path" >&2
                return 1
            }
        done < "$ARTIFACT_AFTER_SNAPSHOT"
    fi

    mv "$ARTIFACT_TMP_MANIFEST" "$manifest_path" || {
        echo "Error: Unable to finalize artifact manifest $manifest_path" >&2
        return 1
    }
    ARTIFACT_TMP_MANIFEST=""
}

cleanup_artifact_tracking() {
    if [ -n "${ARTIFACT_BEFORE_SNAPSHOT:-}" ]; then
        rm -f "$ARTIFACT_BEFORE_SNAPSHOT"
    fi
    if [ -n "${ARTIFACT_AFTER_SNAPSHOT:-}" ]; then
        rm -f "$ARTIFACT_AFTER_SNAPSHOT"
    fi
    if [ -n "${ARTIFACT_TMP_MANIFEST:-}" ]; then
        rm -f "$ARTIFACT_TMP_MANIFEST"
    fi
}

# Check whether a usable genomic FASTA already exists in org_dir for given accession
has_fasta_files() {
    local org_dir="$1"
    local accession="$2"
    [ -n "$(find_first_fasta "$org_dir" "$accession")" ]
}

# Return 0 if f is a non-empty file whose first byte is '>'.
is_valid_fasta() {
    local f="$1"
    [ -s "$f" ] && [[ "$(head -c 1 "$f")" == ">" ]]
}

# Find the first genomic FASTA file in org_dir for given accession.
# Only accepts files matching the accession-prefixed pattern from config and
# excludes non-genomic companions (cds_from_genomic.fna, rna.fna, protein.faa).
# Accepting a non-genomic FASTA as the reference would silently break the
# downstream CDS filter because its seq names wouldn't match the GFF.
find_first_fasta() {
    local org_dir="$1"
    local accession="$2"
    local pattern="${FASTA_PATTERN_TMPL//\{org_id\}/$accession}"
    local f
    while IFS= read -r f; do
        [ -z "$f" ] && continue
        case "$(basename "$f")" in
            cds_from_genomic.fna|rna.fna|protein.faa)
                continue
                ;;
            *_cds_from_genomic.fna|*_rna.fna|*_protein.faa)
                continue
                ;;
        esac
        echo "$f"
        return
    done < <(compgen -G "$org_dir/$pattern" 2>/dev/null || true)
}

# Fetch NCBI Datasets genome summary JSON, parse organism names, save to dest_dir.
# Outputs: final_dir_name|dest_json_path|raw_organism_name|ncbi_full_name|tax_id
get_org_info_from_ncbi() {
    local accession="$1"
    local override_name="$2"
    local dest_dir="$3"

    if ! command -v jq >/dev/null 2>&1; then
        echo "Error: 'jq' is required but not found in PATH." >&2
        return 1
    fi

    local tmp_json
    tmp_json=$(mktemp) || {
        echo "Error: Unable to create temporary file for $accession summary." >&2
        return 1
    }

    if ! curl_with_retry -sS --fail "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$accession/dataset_report" > "$tmp_json"; then
        echo "Error: Failed to fetch genome summary for $accession" >&2
        rm -f "$tmp_json"
        return 1
    fi

    # Parse organism name
    local raw_base_name
    raw_base_name=$(jq -r 'try .reports[0].organism.organism_name catch ""' "$tmp_json")
    if [ -z "$raw_base_name" ] || [ "$raw_base_name" = "null" ]; then
        echo "Warning: organism_name not found in summary for $accession; using accession as name." >&2
        raw_base_name="$accession"
    fi
    local base_name
    base_name=$(sanitize_for_path "$raw_base_name")

    # Parse infraspecific names (strain, cultivar, etc.)
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

    # Parse tax_id
    local tax_id
    tax_id=$(jq -r '(.reports[0].organism.tax_id // .reports[0].biosample.description.organism.tax_id // empty) | tostring' "$tmp_json" 2>/dev/null)

    # Save summary JSON to dest_dir
    if ! mkdir -p "$dest_dir"; then
        echo "Error: Unable to create directory $dest_dir" >&2
        rm -f "$tmp_json"
        return 1
    fi
    local dest_json="$dest_dir/${accession}.json"
    if ! mv "$tmp_json" "$dest_json" 2>/dev/null; then
        cp "$tmp_json" "$dest_json" || {
            echo "Error: Unable to save summary JSON for $accession" >&2
            rm -f "$tmp_json"
            return 1
        }
        rm -f "$tmp_json"
    fi

    echo "$final_dir_name|$dest_json|$raw_base_name|$calculated_full_name|${tax_id:-}"
}

# Unzip ncbi_dataset.zip in org_dir and move genome files up
process_genome_zip() (
    # Resolve to absolute path so it stays valid after cd
    local org_dir
    org_dir="$(cd "$1" && pwd)" || {
        echo "Error: Cannot resolve directory $1" >&2
        return 1
    }
    local accession="$2"

    cd "$org_dir" || {
        echo "Error: Cannot cd to $org_dir" >&2
        return 1
    }

    if [ ! -s "ncbi_dataset.zip" ]; then
        echo "Error: ncbi_dataset.zip is missing or empty in $org_dir" >&2
        return 1
    fi

    if ! unzip ncbi_dataset.zip; then
        echo "Error: Failed to unzip ncbi_dataset.zip in $org_dir" >&2
        return 1
    fi

    cd ncbi_dataset/data || {
        echo "Error: Cannot cd to ncbi_dataset/data in $org_dir" >&2
        return 1
    }

    # Move top-level files (not directories)
    if ! mv $(ls -p | grep -v /) "$org_dir" 2>/dev/null; then
        true  # top-level files may not exist; continue
    fi

    if [ -d "$accession" ]; then
        mv "$accession"/* "$org_dir" || {
            echo "Error: Failed to move $accession/* to $org_dir" >&2
            return 1
        }
    fi

    cd "$org_dir" || return 1
    rm -rf ncbi_dataset || echo "Warning: Could not remove ncbi_dataset directory" >&2
    echo "Successfully processed genome data for accession $accession"
)

# --- argument parsing ---

ACCESSION=""
OUT_DIR="."
OVERRIDE_NAME=""
NO_GENOME=0
NO_TAXONOMY=0
TAX_OUT_OVERRIDE=""
ARTIFACT_MANIFEST=""
ARTIFACT_BEFORE_SNAPSHOT=""
ARTIFACT_AFTER_SNAPSHOT=""
ARTIFACT_TMP_MANIFEST=""
FORCE=0

trap cleanup_artifact_tracking EXIT

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir)
            OUT_DIR="${2:?--out-dir requires a path}"
            shift 2
            ;;
        --out-dir=*)
            OUT_DIR="${1#*=}"
            [ -z "$OUT_DIR" ] && { echo "Error: --out-dir requires a path" >&2; exit 1; }
            shift
            ;;
        --name)
            OVERRIDE_NAME="${2:?--name requires a value}"
            shift 2
            ;;
        --name=*)
            OVERRIDE_NAME="${1#*=}"
            shift
            ;;
        --no-genome)
            NO_GENOME=1
            shift
            ;;
        --no-taxonomy)
            NO_TAXONOMY=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --tax-out)
            TAX_OUT_OVERRIDE="${2:?--tax-out requires a path}"
            shift 2
            ;;
        --tax-out=*)
            TAX_OUT_OVERRIDE="${1#*=}"
            shift
            ;;
        --artifact-manifest)
            ARTIFACT_MANIFEST="${2:?--artifact-manifest requires a path}"
            shift 2
            ;;
        --artifact-manifest=*)
            ARTIFACT_MANIFEST="${1#*=}"
            if [ -z "$ARTIFACT_MANIFEST" ]; then
                echo "Error: --artifact-manifest requires a path" >&2
                exit 1
            fi
            shift
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)
            if [ -z "$ACCESSION" ]; then
                ACCESSION="$1"
            else
                echo "Error: Unexpected argument $1" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$ACCESSION" ]; then
    echo "Usage: $(basename "$0") <ACCESSION> [--out-dir DIR] [--name NAME] [--no-genome] [--no-taxonomy] [--tax-out FILE] [--artifact-manifest FILE] [--force]" >&2
    exit 1
fi

# Load fasta pattern from config (fall back to default)
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
config_file="$ROOT_DIR/config/dwl_config.yaml"
if command -v yq >/dev/null 2>&1 && [ -f "$config_file" ]; then
    FASTA_PATTERN_TMPL=$(yq eval '.patterns.fasta' "$config_file" 2>/dev/null)
fi
[ -z "$FASTA_PATTERN_TMPL" ] && FASTA_PATTERN_TMPL='{org_id}*.fna'

# Build download query params from config annotation types
DWL_QUERY_PARAMS=""
if command -v yq >/dev/null 2>&1 && [ -f "$config_file" ]; then
    while IFS= read -r annot_type; do
        [ -n "$annot_type" ] && DWL_QUERY_PARAMS="${DWL_QUERY_PARAMS}&include_annotation_type=${annot_type}"
    done < <(yq eval '.download.includes[]' "$config_file" 2>/dev/null)
fi
[ -z "$DWL_QUERY_PARAMS" ] && DWL_QUERY_PARAMS="include_annotation_type=GENOME_FASTA"
DWL_QUERY_PARAMS="${DWL_QUERY_PARAMS#&}"

# --- step 1: fetch summary JSON and resolve organism name ---

mkdir -p "$OUT_DIR"
# Temporary dest_dir placeholder; real dir determined after name resolution
tmp_dest="$OUT_DIR/__tmp_dwl_${ACCESSION}"
mkdir -p "$tmp_dest"

ncbi_info=$(get_org_info_from_ncbi "$ACCESSION" "$OVERRIDE_NAME" "$tmp_dest") || exit 1

IFS='|' read -r dir_name summary_json_path raw_organism_name ncbi_full_name tax_id <<< "$ncbi_info"

# Move summary JSON to final org dir
org_dir="$OUT_DIR/$dir_name"
org_dir_preexisting=0
if [ -n "$ARTIFACT_MANIFEST" ] && { [ -e "$org_dir" ] || [ -L "$org_dir" ]; }; then
    org_dir_preexisting=1
    ARTIFACT_BEFORE_SNAPSHOT=$(mktemp) || {
        echo "Error: Unable to create artifact snapshot." >&2
        exit 1
    }
    snapshot_relative_paths "$org_dir" "$ARTIFACT_BEFORE_SNAPSHOT" || {
        echo "Error: Unable to snapshot existing artifacts in $org_dir" >&2
        exit 1
    }
fi
mkdir -p "$org_dir"
final_summary_json="$org_dir/${ACCESSION}.json"
if [ "$summary_json_path" != "$final_summary_json" ]; then
    mv "$summary_json_path" "$final_summary_json" 2>/dev/null || cp "$summary_json_path" "$final_summary_json"
fi
rmdir "$tmp_dest" 2>/dev/null || true
summary_json_path="$final_summary_json"

echo "Organism: $ncbi_full_name (dir: $dir_name)" >&2

# --- step 2: genome download ---

fasta_path=""
if [ "$NO_GENOME" -eq 0 ]; then
    _need_download=1
    if [ "$FORCE" -eq 0 ] && has_fasta_files "$org_dir" "$ACCESSION"; then
        _existing=$(find_first_fasta "$org_dir" "$ACCESSION")
        if is_valid_fasta "$_existing"; then
            echo "Genome files already exist for $dir_name; skipping download." >&2
            _need_download=0
        else
            echo "Existing FASTA failed validation; re-downloading $dir_name ($ACCESSION)..." >&2
        fi
    fi
    if [ "$_need_download" -eq 1 ]; then
        echo "Downloading genome for $dir_name ($ACCESSION)..." >&2
        (
            cd "$org_dir"
            if ! curl_with_retry -sS --fail -o ncbi_dataset.zip \
                "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$ACCESSION/download?${DWL_QUERY_PARAMS}"; then
                echo "Error: Failed to download genome for $ACCESSION" >&2
                exit 1
            fi
            if [ ! -s "ncbi_dataset.zip" ]; then
                echo "Error: ncbi_dataset.zip is empty or missing after download" >&2
                exit 1
            fi
        ) || exit 1
        process_genome_zip "$org_dir" "$ACCESSION" >&2 || exit 1
    fi
    fasta_path=$(find_first_fasta "$org_dir" "$ACCESSION")
    if [ -z "$fasta_path" ]; then
        echo "Warning: No FASTA file found in $org_dir after download." >&2
    fi
fi

# --- step 3: taxonomy download ---

tax_json_path=""
if [ "$NO_TAXONOMY" -eq 0 ]; then
    if [ -n "$TAX_OUT_OVERRIDE" ]; then
        tax_json_path="$TAX_OUT_OVERRIDE"
    else
        tax_json_path="$org_dir/taxonomy_${ACCESSION}.json"
    fi

    if [ -s "$tax_json_path" ]; then
        echo "Taxonomy JSON already exists: $tax_json_path" >&2
    elif [ -z "$tax_id" ] || [ "$tax_id" = "null" ]; then
        echo "Warning: tax_id not found in summary for $ACCESSION; skipping taxonomy download." >&2
        tax_json_path=""
    else
        echo "Downloading taxonomy for tax_id=$tax_id..." >&2
        tmp_tax=$(mktemp) || { echo "Warning: Cannot create tmpfile for taxonomy." >&2; tax_json_path=""; }
        if [ -n "$tmp_tax" ]; then
            if curl_with_retry -sS --fail "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/$tax_id/dataset_report" > "$tmp_tax" && [ -s "$tmp_tax" ]; then
                mv "$tmp_tax" "$tax_json_path" 2>/dev/null || {
                    cp "$tmp_tax" "$tax_json_path"
                    rm -f "$tmp_tax"
                }
                echo "Downloaded taxonomy JSON: $tax_json_path" >&2
            else
                echo "Warning: taxonomy download failed or returned empty for tax_id=$tax_id" >&2
                rm -f "$tmp_tax"
                tax_json_path=""
            fi
        fi
    fi
fi

# --- optional artifact manifest ---

if [ -n "$ARTIFACT_MANIFEST" ]; then
    write_artifact_manifest \
        "$ARTIFACT_MANIFEST" "$org_dir" "$dir_name" "$org_dir_preexisting" \
        "$ARTIFACT_BEFORE_SNAPSHOT" || exit 1
fi

# --- stdout: pipe-separated result ---

echo "${dir_name}|${fasta_path}|${summary_json_path}|${raw_organism_name}|${ncbi_full_name}|${tax_json_path}"
