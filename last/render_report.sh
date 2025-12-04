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
OUTPUT_NAME=""
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

DATASET_NAMES="$(
python - "$INPUT_JSON_ABS" <<'PY'
import json, sys
path = sys.argv[1]
with open(path, "r", encoding="utf-8") as handle:
    data = json.load(handle)
datasets = data.get("datasets") or []
if not datasets:
    print("Error: summary JSON contains no datasets.", file=sys.stderr)
    sys.exit(1)
required_keys = ("species1", "species2", "species3", "idt_12", "idt_13", "idt_23")
for entry in datasets:
    name = entry.get("dataset") or "unknown dataset"
    missing = [key for key in required_keys if key not in entry]
    if missing:
        print(
            f"Error: {name} is missing required fields from collect_run_summary.py: {', '.join(missing)}",
            file=sys.stderr,
        )
        sys.exit(1)
print(",".join(entry.get("dataset", "unknown") for entry in datasets))
PY
)"
echo "Rendering datasets: $DATASET_NAMES"

if [[ -n "$OUTPUT_NAME" ]]; then
    OUTPUT_PATH="$(python - "$OUTPUT_NAME" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"
else
    OUTPUT_PATH="$(python - "$INPUT_JSON_ABS" <<'PY'
import os, sys
json_path = os.path.abspath(sys.argv[1])
base = os.path.splitext(os.path.basename(json_path))[0] or "report"
default_docx = base + ".docx"
print(os.path.join(os.path.dirname(json_path), default_docx))
PY
)"
fi

PREVIEW_DIR="$(python - "$INPUT_JSON_ABS" <<'PY'
import os, sys
json_dir = os.path.dirname(os.path.abspath(sys.argv[1]))
preview_dir = os.path.join(json_dir, "tmp")
os.makedirs(preview_dir, exist_ok=True)
print(preview_dir)
PY
)"

IMAGES_BASE_DIR="$SCRIPT_DIR/images"
mkdir -p "$IMAGES_BASE_DIR"
IMAGES_DIR="$(mktemp -d "$IMAGES_BASE_DIR/run.XXXXXX")"

cleanup_images() {
    rm -rf "$IMAGES_DIR"
}
trap cleanup_images EXIT

Rscript --vanilla - "$RMD_PATH" "$INPUT_JSON_ABS" "$OUTPUT_PATH" "$OUTPUT_FORMAT" "$PREVIEW_DIR" "$IMAGES_DIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
rmd_path <- args[[1]]
input_json <- args[[2]]
output_file <- args[[3]]
output_format <- args[[4]]
preview_dir <- args[[5]]
images_dir <- args[[6]]

adjust_docx_image_widths <- function(docx_path, target_cm = 16.5) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    message("Skipping width normalization; 'xml2' package not available.")
    return(invisible(FALSE))
  }

  target_emu <- as.integer(round(target_cm / 2.54 * 914400))
  work_dir <- tempfile("docx_fix_")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  utils::unzip(docx_path, exdir = work_dir)
  doc_xml_path <- file.path(work_dir, "word", "document.xml")
  if (!file.exists(doc_xml_path)) {
    return(invisible(FALSE))
  }

  doc <- xml2::read_xml(doc_xml_path)
  ns <- xml2::xml_ns(doc)
  extent_nodes <- xml2::xml_find_all(doc, ".//wp:extent", ns)
  if (length(extent_nodes) == 0) {
    return(invisible(FALSE))
  }

  changed <- FALSE
  for (node in extent_nodes) {
    parent_inline <- xml2::xml_parent(node)
    pic_descr <- xml2::xml_find_first(parent_inline, ".//pic:pic/pic:nvPicPr/pic:cNvPr", ns)
    descr <- if (inherits(pic_descr, "xml_missing")) NA_character_ else xml2::xml_attr(pic_descr, "descr")
    if (is.na(descr) || !startsWith(descr, "images/")) {
      next
    }

    cx_old <- suppressWarnings(as.numeric(xml2::xml_attr(node, "cx")))
    cy_old <- suppressWarnings(as.numeric(xml2::xml_attr(node, "cy")))
    if (!is.finite(cx_old) || cx_old <= 0) {
      next
    }

    scale <- target_emu / cx_old
    xml2::xml_set_attr(node, "cx", as.character(target_emu))
    if (is.finite(cy_old) && cy_old > 0) {
      xml2::xml_set_attr(node, "cy", as.character(round(cy_old * scale)))
    }
    changed <- TRUE
  }

  if (!changed) {
    return(invisible(FALSE))
  }

  xml2::write_xml(doc, doc_xml_path)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(work_dir)
  zip_path <- normalizePath(docx_path, mustWork = FALSE)
  if (file.exists(zip_path)) {
    file.remove(zip_path)
  }
  files_to_zip <- list.files(".", recursive = TRUE, all.files = TRUE)
  zip_status <- system2("zip", c("-r9Xq", zip_path, files_to_zip), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(zip_status, "status")) && attr(zip_status, "status") != 0) {
    stop("Failed to repack DOCX: ", paste(zip_status, collapse = "\n"))
  }
  invisible(TRUE)
}

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

if (!dir.exists(images_dir)) {
  dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)
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
    images_dir = images_dir,
    pdftools_available = pdt_available
  ),
  output_file = output_file,
  output_format = output_format
)

format_str <- tolower(paste(output_format, collapse = ""))
if (nzchar(format_str) && grepl("word_document", format_str, fixed = TRUE)) {
  try(adjust_docx_image_widths(output_file), silent = TRUE)
}
RSCRIPT

echo "Rendered report: $OUTPUT_PATH"


