#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RMD_PATH="$SCRIPT_DIR/report_template.Rmd"
CONTRACT_PY="$SCRIPT_DIR/report_contract.py"
RENDER_R="$SCRIPT_DIR/render_report.R"

usage() {
    cat <<EOF
Usage: $(basename "$0") -j <input_json> [-o output_file] [-f output_format]

Options:
  -j, --json       Path to the JSON summary produced by collect_run_summary.py (required).
  -o, --output     Output filename for the rendered document.
                   Defaults to the JSON basename with a format-appropriate extension.
  -f, --format     rmarkdown output format (default: word_document).
                   Examples: word_document, pdf_document, html_document.
  -h, --help       Show this help message and exit.

Examples:
  $(basename "$0") -j results/fungi/fungi_summary.json
  $(basename "$0") -j results/fungi/fungi_summary.json -o fungi_summary.pdf -f pdf_document
EOF
}

INPUT_JSON=""
OUTPUT_NAME=""
OUTPUT_FORMAT="word_document"

require_option_value() {
  local option="$1"
  local value="${2:-}"
  if [[ -z "$value" || "$value" == -* ]]; then
    echo "Error: ${option} requires a value." >&2
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--json)
          require_option_value "$1" "${2:-}"
          INPUT_JSON="$2"
            shift 2
            ;;
        --json=*)
          INPUT_JSON="${1#*=}"
          require_option_value "--json" "$INPUT_JSON"
          shift
          ;;
        -o|--output)
          require_option_value "$1" "${2:-}"
          OUTPUT_NAME="$2"
            shift 2
            ;;
        --output=*)
          OUTPUT_NAME="${1#*=}"
          require_option_value "--output" "$OUTPUT_NAME"
          shift
          ;;
        -f|--format)
          require_option_value "$1" "${2:-}"
          OUTPUT_FORMAT="$2"
            shift 2
            ;;
        --format=*)
          OUTPUT_FORMAT="${1#*=}"
          require_option_value "--format" "$OUTPUT_FORMAT"
          shift
          ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_JSON" ]]; then
    echo "Error: --json is required." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$INPUT_JSON" ]]; then
    echo "Error: JSON file not found: $INPUT_JSON" >&2
    exit 1
fi

if [[ ! -f "$RMD_PATH" ]]; then
    echo "Error: R Markdown template not found at ${RMD_PATH}." >&2
    exit 1
fi

if [[ ! -f "$CONTRACT_PY" ]]; then
  echo "Error: report_contract.py not found at ${CONTRACT_PY}." >&2
  exit 1
fi

if [[ ! -f "$RENDER_R" ]]; then
  echo "Error: render_report.R not found at ${RENDER_R}." >&2
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "Error: python3 is not available in PATH." >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "Error: Rscript is not available in PATH." >&2
    exit 1
fi

python3 "$CONTRACT_PY" validate-format "$OUTPUT_FORMAT"

INPUT_JSON_ABS="$(realpath "$INPUT_JSON")"
DATASET_NAMES="$(python3 "$CONTRACT_PY" validate "$INPUT_JSON_ABS")"
echo "Rendering datasets: $DATASET_NAMES"

OUTPUT_PATH_ARGS=("report-path" "$INPUT_JSON_ABS" "$OUTPUT_FORMAT")
if [[ -n "$OUTPUT_NAME" ]]; then
  OUTPUT_PATH_ARGS+=("--output" "$OUTPUT_NAME")
fi
OUTPUT_PATH="$(python3 "$CONTRACT_PY" "${OUTPUT_PATH_ARGS[@]}")"

PREVIEW_DIR="$(python3 "$CONTRACT_PY" preview-dir "$INPUT_JSON_ABS")"

IMAGES_BASE_DIR="$SCRIPT_DIR/images"
mkdir -p "$IMAGES_BASE_DIR"
IMAGES_DIR="$(mktemp -d "$IMAGES_BASE_DIR/run.XXXXXX")"

cleanup_images() {
    rm -rf "$IMAGES_DIR"
}
trap cleanup_images EXIT

PDFTOOLS_AVAIL="FALSE"
if Rscript --vanilla \
    -e "quit(status = if (requireNamespace('pdftools', quietly=TRUE)) 0L else 1L)" \
    >/dev/null 2>&1; then
    PDFTOOLS_AVAIL="TRUE"
fi

Rscript --vanilla "$RENDER_R" \
    "$RMD_PATH" \
    "$INPUT_JSON_ABS" \
    "$OUTPUT_PATH" \
    "$OUTPUT_FORMAT" \
    "$PREVIEW_DIR" \
    "$IMAGES_DIR" \
    "$PDFTOOLS_AVAIL"

echo "Rendered report: $OUTPUT_PATH"


