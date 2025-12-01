#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RMD_PATH="$SCRIPT_DIR/report_template.Rmd"

usage() {
    cat <<EOF
Usage: $(basename "$0") -j <input_json> [-o output_file] [-f output_format]

Options:
  -j, --json       Path to the JSON summary produced by collect_run_summary.py (required).
  -o, --output     Output filename for the rendered document (default: fungi_summary.docx).
  -f, --format     rmarkdown output format (default: word_document).
                   Examples: word_document, pdf_document, html_document.
  -h, --help       Show this help message and exit.

Examples:
  $(basename "$0") -j /home/mrk/sbst/data/fungi/fungi_summary.json
  $(basename "$0") -j /home/mrk/sbst/data/fungi/fungi_summary.json -o fungi_summary.pdf -f pdf_document
EOF
}

INPUT_JSON=""
OUTPUT_NAME="fungi_summary.docx"
OUTPUT_FORMAT="word_document"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--json)
            INPUT_JSON="${2:-}"
            shift 2
            ;;
        -o|--output)
            OUTPUT_NAME="${2:-}"
            shift 2
            ;;
        -f|--format)
            OUTPUT_FORMAT="${2:-}"
            shift 2
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

if ! command -v Rscript >/dev/null 2>&1; then
    echo "Error: Rscript is not available in PATH." >&2
    exit 1
fi

INPUT_JSON_ABS="$(python - "$INPUT_JSON" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"

OUTPUT_PATH="$(python - "$OUTPUT_NAME" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"

PREVIEW_DIR="$(mktemp -d)"
cleanup() {
    rm -rf "$PREVIEW_DIR"
}
trap cleanup EXIT

Rscript --vanilla - "$RMD_PATH" "$INPUT_JSON_ABS" "$OUTPUT_PATH" "$OUTPUT_FORMAT" "$PREVIEW_DIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
rmd_path <- args[[1]]
input_json <- args[[2]]
output_file <- args[[3]]
output_format <- args[[4]]
preview_dir <- args[[5]]

if (!nzchar(output_format)) {
  output_format <- NULL
}

if (!file.exists(rmd_path)) {
  stop(sprintf("Template not found: %s", rmd_path))
}

if (!file.exists(input_json)) {
  stop(sprintf("Input JSON not found: %s", input_json))
}

if (!dir.exists(preview_dir)) {
  dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
}

packages <- c("rmarkdown")
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf("Missing required R packages: %s", paste(missing_pkgs, collapse = ", ")))
}

if (!requireNamespace("pdftools", quietly = TRUE)) {
  message("Attempting to install 'pdftools' package for PDF preview support...")
  utils::install.packages("pdftools", repos = "https://cloud.r-project.org")
}

pdt_available <- requireNamespace("pdftools", quietly = TRUE)

rmarkdown::render(
  input = rmd_path,
  params = list(
    input_json = input_json,
    preview_dir = preview_dir,
    pdftools_available = pdt_available
  ),
  output_file = output_file,
  output_format = output_format
)
RSCRIPT

echo "Rendered report: $OUTPUT_PATH"


