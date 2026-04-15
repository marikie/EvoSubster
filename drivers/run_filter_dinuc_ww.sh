#! /bin/bash
# Apply filter_top_n_context_dinuc.R for three W-W complement-pair queries
# against every chi2_context_dinuc output file in results/fungi/chi2_dinuc_glom/.
#
# Queries:
#   TT  TT>AA
#   TA  TA>AT
#   AT  AT>TA
#
# Each file becomes one === section === in a temporary combined input.
# All three query results are written sequentially to a single output file.
#
# Usage: bash drivers/run_filter_dinuc_ww.sh [--out FILE]
#   --out FILE  Output path (default: results/fungi/chi2_dinuc_glom/filter_ww.txt)

set -euo pipefail

DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$DRIVER_DIR/.." && pwd)"
R_SCRIPT="$ROOT_DIR/src/statistics/filter_top_n_context_dinuc.R"
INPUT_DIR="$ROOT_DIR/results/fungi/chi2_dinuc_glom"
OUT_FILE="$INPUT_DIR/filter_ww.txt"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT_FILE="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Build combined input: each *_ncds.txt and *_full.txt becomes === name === section
COMBINED=$(mktemp)
trap 'rm -f "$COMBINED"' EXIT

for f in "$INPUT_DIR"/*_ncds.txt "$INPUT_DIR"/*_full.txt; do
  [[ -f "$f" ]] || continue
  label=$(basename "${f%.txt}")
  echo "=== $label ===" >> "$COMBINED"
  cat "$f"              >> "$COMBINED"
done

# Queries: "doublet context"
QUERIES=(
  "TT TT>AA"
  "TA TA>AT"
  "AT AT>TA"
)

> "$OUT_FILE"

for query in "${QUERIES[@]}"; do
  read -r doublet ctx <<< "$query"
  {
    echo "========================================"
    echo "# Filter: $doublet  $ctx"
    echo "========================================"
    Rscript "$R_SCRIPT" "$COMBINED" "$doublet" "$ctx"
  } >> "$OUT_FILE"
done

echo "Done. Output: $OUT_FILE"
