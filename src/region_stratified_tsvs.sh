#!/bin/bash
# region_stratified_tsvs.sh
#
# Given a three-way joined MAF and a region annotation (BED/GFF) on the TOP
# (outgroup) sequence, produce substitution-count TSVs separately for the
# region class (inside, e.g. TE) and its complement (outside, e.g. non-TE),
# for both ingroups and for both single-base (96 SBS) and double-base spectra.
#
# This is the region-stratified generalisation of generate_tsv_files.sh:
# it wires the new src/align/maf-cut-region.py splitter into the existing,
# unchanged counters src/count/single_sbst_2TSVs.py and disbst_2TSVs.py.
#
# Usage:
#   src/region_stratified_tsvs.sh <joined.maf> <regions.bed|gff> <out_dir> \
#       [--label NAME] [--type GFF_TYPE] [--shrink]
#
#   <joined.maf>   org1_org2_org3_*.maf (top sequence = outgroup)
#   <regions>      BED or GFF of regions on the outgroup assembly
#   <out_dir>      output directory (created if absent)
#   --label NAME   class label used in filenames (default: region)
#   --type T       GFF feature type filter (e.g. CDS); ignored for BED
#   --shrink       shrink regions 1 bp/side (maf-cut-cds-uglier.py compatibility)
#
# Outputs under <out_dir>/:
#   <label>/singlenuc/{org2,org3}.tsv        non<label>/singlenuc/{org2,org3}.tsv
#   <label>/dinuc/{org2,org3}_dinuc.tsv      non<label>/dinuc/{org2,org3}_dinuc.tsv
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CUTTER="$SCRIPT_DIR/align/maf-cut-region.py"
SBS="$SCRIPT_DIR/count/single_sbst_2TSVs.py"
DBS="$SCRIPT_DIR/count/disbst_2TSVs.py"

if [ $# -lt 3 ]; then
    sed -n '15,33p' "$0"; exit 2
fi

MAF="$1"; REGIONS="$2"; OUT="$3"; shift 3
LABEL="region"; GFFTYPE=""; SHRINK=""
while [ $# -gt 0 ]; do
    case "$1" in
        --label) LABEL="$2"; shift 2 ;;
        --type)  GFFTYPE="$2"; shift 2 ;;
        --shrink) SHRINK="--shrink"; shift ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

[ -f "$MAF" ]     || { echo "no such MAF: $MAF" >&2; exit 1; }
[ -f "$REGIONS" ] || { echo "no such regions file: $REGIONS" >&2; exit 1; }

# derive org2/org3 labels from filename: org1_org2_org3_*.maf
base="$(basename "$MAF")"; stem="${base%.maf}"
IFS='_' read -r -a parts <<< "$stem"
if [ "${#parts[@]}" -lt 3 ]; then
    echo "MAF name must be org1_org2_org3_*.maf, got: $base" >&2; exit 1
fi
ORG2="${parts[1]}"; ORG3="${parts[2]}"

typeflag=()
[ -n "$GFFTYPE" ] && typeflag=(--type "$GFFTYPE")

# inside (region) = $LABEL ; outside (complement) = non$LABEL
run_one() {  # $1 keep(inside|outside)  $2 outlabel
    local keep="$1" outlabel="$2"
    local tmp; tmp="$(mktemp -d)"
    local cut="$tmp/${outlabel}.maf"
    python3 "$CUTTER" "$REGIONS" "$MAF" --keep "$keep" "${typeflag[@]}" $SHRINK > "$cut"

    mkdir -p "$OUT/$outlabel/singlenuc" "$OUT/$outlabel/dinuc"
    python3 "$SBS" "$cut" \
        -o2 "$OUT/$outlabel/singlenuc/${ORG2}.tsv" \
        -o3 "$OUT/$outlabel/singlenuc/${ORG3}.tsv"
    python3 "$DBS" "$cut" \
        -o2 "$OUT/$outlabel/dinuc/${ORG2}_dinuc.tsv" \
        -o3 "$OUT/$outlabel/dinuc/${ORG3}_dinuc.tsv"
    rm -rf "$tmp"
    echo "  wrote $OUT/$outlabel/{singlenuc,dinuc}/"
}

echo "Region-stratified counting for $base"
echo "  regions: $REGIONS (label=$LABEL)"
run_one inside  "$LABEL"
run_one outside "non$LABEL"
echo "done."
