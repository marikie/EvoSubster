#!/bin/bash
# Run chi2_context_dinuc.R on all Glomeromycetes trios (trio_glom.md).
#
# Usage: bash run_chi2_dinuc_glom.sh [--out-dir DIR]
#   --out-dir DIR  Directory for output files (default: results/fungi/chi2_dinuc_glom)
#
# For Non-coding:Yes trios, uses _maflinked_dinuc_ncds.tsv.
# For Non-coding:No  trios, uses _maflinked_dinuc.tsv.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RESULTS_DIR="$REPO_ROOT/results/fungi"
R_SCRIPT="$SCRIPT_DIR/chi2_context_dinuc.R"

# Default output directory
OUT_DIR="$RESULTS_DIR/chi2_dinuc_glom"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"

# Trio name -> non-coding flag (1 = use _ncds, 0 = full genome)
declare -A NON_CODING
NON_CODING["Denhet1_Gigros2_Gigmar3"]=1
NON_CODING["Fungeo1_Funcal2_Funmos3"]=1
NON_CODING["Parbra1_Parocc2_Parocc3"]=1
NON_CODING["Rhicla1_Rhiirr2_Rhipro3"]=1
NON_CODING["Rhipro1_Rhiirr2_Rhiirr3"]=0

TRIOS=(
  Denhet1_Gigros2_Gigmar3
  Fungeo1_Funcal2_Funmos3
  Parbra1_Parocc2_Parocc3
  Rhicla1_Rhiirr2_Rhipro3
  Rhipro1_Rhiirr2_Rhiirr3
)

for trio in "${TRIOS[@]}"; do
  trio_dir="$RESULTS_DIR/$trio"
  date_dir=$(ls "$trio_dir" | head -1)
  run_dir="$trio_dir/$date_dir"
  ncds="${NON_CODING[$trio]}"

  if [[ "$ncds" -eq 1 ]]; then
    pattern="*_maflinked_dinuc_ncds.tsv"
    suffix="ncds"
  else
    pattern="*_maflinked_dinuc.tsv"
    suffix="full"
  fi

  mapfile -t tsvs < <(find "$run_dir" -maxdepth 1 -name "$pattern" | sort)

  if [[ ${#tsvs[@]} -eq 0 ]]; then
    echo "WARNING: No TSV files matching '$pattern' in $run_dir" >&2
    continue
  fi

  for tsv in "${tsvs[@]}"; do
    filename=$(basename "$tsv")
    # Extract species abbreviation from filename: GCA_ACC_SpeciesAbbr_DATE_...
    species_abbr=$(echo "$filename" | awk -F'_' '{print $3}')
    label="${trio} / ${species_abbr} (${suffix})"
    out_file="$OUT_DIR/${trio}_${species_abbr}_${suffix}.txt"

    echo "Running: $trio / $species_abbr [$suffix]"
    Rscript "$R_SCRIPT" "$tsv" "$label" > "$out_file"
    echo "  -> $out_file"
  done
done

echo
echo "Done. Results in: $OUT_DIR"
