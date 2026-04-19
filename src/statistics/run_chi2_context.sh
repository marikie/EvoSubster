#!/bin/bash
# run_chi2_context.sh - Run chi2_context.R on all trios in a lineage directory.
#
# For each trio, finds org2 and org3 TSVs (*_ncds.tsv preferred,
# fallback to *.tsv), runs chi2_context.R, and aggregates results.
#
# Usage: bash run_chi2_context.sh <lineage_dir> [output_file]
#   lineage_dir  : e.g. results/fungi/
#   output_file  : default <lineage_dir>/<lineage_name>_chi2_context.txt

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <lineage_dir> [output_file]" >&2
    exit 1
fi

LINEAGE_DIR="${1%/}"  # strip trailing slash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/chi2_context.R"

if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: chi2_context.R not found at ${R_SCRIPT}" >&2
    exit 1
fi

LINEAGE_NAME="$(basename "$LINEAGE_DIR")"

if [ $# -ge 2 ]; then
    OUTPUT_FILE="$2"
else
    OUTPUT_FILE="${LINEAGE_DIR}/${LINEAGE_NAME}_chi2_context.txt"
fi

# Write header
{
    echo "# Chi2 context dependence -- ${LINEAGE_NAME} ($(date +%Y-%m-%d))"
    echo "# Input: ${LINEAGE_DIR}"
    echo ""
} > "$OUTPUT_FILE"

find_tsv() {
    local date_dir="$1"
    local org_label="$2"
    local files

    local singlenuc_dir="${date_dir}/statistics/${org_label}/singlenuc"

    # Prefer _ncds; use array glob so nullglob yields empty array on no match
    files=("${singlenuc_dir}"/*_ncds.tsv)
    if [ "${#files[@]}" -gt 0 ] && [ -f "${files[0]}" ]; then
        echo "${files[0]}"
        return
    fi

    # Fallback to non-ncds
    files=("${singlenuc_dir}"/*.tsv)
    if [ "${#files[@]}" -gt 0 ] && [ -f "${files[0]}" ]; then
        echo "${files[0]}"
        return
    fi

    echo ""
}

run_for_org() {
    local trio_name="$1"
    local date_dir="$2"
    local org_label="$3"

    local tsv
    tsv=$(find_tsv "$date_dir" "$org_label")

    if [ -z "$tsv" ]; then
        echo "WARNING: no TSV found for ${trio_name}/${org_label} in ${date_dir}" >&2
        return
    fi

    local label="${trio_name}/${org_label}"
    echo "=== ${label} ===" >> "$OUTPUT_FILE"
    Rscript "$R_SCRIPT" "$tsv" "$label" >> "$OUTPUT_FILE" 2>&1
    echo "" >> "$OUTPUT_FILE"
}

shopt -s nullglob

for trio_dir in "${LINEAGE_DIR}"/*/; do
    [ -d "$trio_dir" ] || continue
    trio_name="$(basename "$trio_dir")"

    # Extract org2 and org3 from trio directory name (format: Org1_Org2_Org3)
    org2="$(echo "$trio_name" | cut -d'_' -f2)"
    org3="$(echo "$trio_name" | cut -d'_' -f3)"

    if [ -z "$org2" ] || [ -z "$org3" ]; then
        echo "WARNING: cannot parse trio name '${trio_name}', skipping" >&2
        continue
    fi

    # Use latest date directory
    date_dir=$(ls -d "${trio_dir}"*/ 2>/dev/null | sort | tail -1)
    if [ -z "$date_dir" ]; then
        echo "WARNING: no date directory found under ${trio_dir}, skipping" >&2
        continue
    fi

    run_for_org "$trio_name" "$date_dir" "$org2"
    run_for_org "$trio_name" "$date_dir" "$org3"
done

echo "Done. Output: ${OUTPUT_FILE}"
