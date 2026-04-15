#! /bin/bash
# Run chi2_context_dinuc.R on all Glomeromycetes trios (trio_glom.md).
#
# For Non-coding:Yes trios, processes both _maflinked_dinuc_ncds.tsv
# and _maflinked_dinuc.tsv. For Non-coding:No trios, processes only
# _maflinked_dinuc.tsv.
#
# Usage: bash drivers/run_chi2_dinuc_glom.sh [--out-dir DIR]
#   --out-dir DIR  Directory for output files
#                  (default: results/fungi/chi2_dinuc_glom)

set -euo pipefail

DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$DRIVER_DIR/.." && pwd)"
R_SCRIPT="$ROOT_DIR/src/statistics/chi2_context_dinuc.R"
RESULTS_DIR="$ROOT_DIR/results/fungi"
OUT_DIR="$RESULTS_DIR/chi2_dinuc_glom"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"

# Trio name -> non-coding flag (1 = also run _ncds variant, 0 = full only)
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

  # Build list of (pattern, suffix) pairs to process
  if [[ "$ncds" -eq 1 ]]; then
    patterns=("*_maflinked_dinuc_ncds.tsv" "*_maflinked_dinuc.tsv")
    suffixes=("ncds" "full")
  else
    patterns=("*_maflinked_dinuc.tsv")
    suffixes=("full")
  fi

  for i in "${!patterns[@]}"; do
    pattern="${patterns[$i]}"
    suffix="${suffixes[$i]}"

    mapfile -t tsvs < <(find "$run_dir" -maxdepth 1 -name "$pattern" | sort)

    if [[ ${#tsvs[@]} -eq 0 ]]; then
      echo "WARNING: No TSV files matching '$pattern' in $run_dir" >&2
      continue
    fi

    for tsv in "${tsvs[@]}"; do
      filename=$(basename "$tsv")
      # Filename format: GCA_ACC_SpeciesAbbr_DATE_...
      species_abbr=$(echo "$filename" | awk -F'_' '{print $3}')
      label="${trio} / ${species_abbr} (${suffix})"
      out_file="$OUT_DIR/${trio}_${species_abbr}_${suffix}.txt"

      echo "Running: $trio / $species_abbr [$suffix]"
      Rscript "$R_SCRIPT" "$tsv" "$label" > "$out_file"
      echo "  -> $out_file"
    done
  done
done

echo
echo "Done. Results in: $OUT_DIR"
