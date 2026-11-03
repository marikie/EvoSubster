#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COLLECT_PY="$SCRIPT_DIR/collect_run_summary.py"
CONTRACT_PY="$SCRIPT_DIR/report_contract.py"
RENDER_SH="$SCRIPT_DIR/render_report.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") <input_dir> [options]

Run collect_run_summary.py then render_report.sh in sequence.

Arguments:
    input_dir                   Lineage root (results/<lineage>) or a single trio root.

Options:
  --idt-threshold FLOAT       Identity threshold for dataset filtering (default: 80.0).
  --json PATH                 Output path for the JSON summary.
                              Defaults to <input_dir>/<name>_summary.json.
  --use-filtered              Use the filtered JSON (_filtered.json) for report rendering.
  -o, --report-output PATH    Output filename for the rendered report.
  -f, --format FORMAT         rmarkdown output format (default: word_document).
                              Examples: word_document, html_document.
  --collect-only              Run only collect_run_summary.py; skip report rendering.
  -h, --help                  Show this help message and exit.

Examples:
    $(basename "$0") results/fungi
    $(basename "$0") results/fungi --use-filtered
    $(basename "$0") results/fungi --idt-threshold 85 --use-filtered
    $(basename "$0") results/Oikdio1_Oikalb2_Oikvan3 -f html_document -o /tmp/report.html
    $(basename "$0") results/fungi --collect-only
EOF
}

IDT_THRESHOLD=""
JSON_OUTPUT=""
USE_FILTERED=0
REPORT_OUTPUT=""
OUTPUT_FORMAT="word_document"
COLLECT_ONLY=0
POSITIONAL_ARGS=()

require_option_value() {
    local option="$1"
    local value="${2:-}"
    if [[ -z "$value" || ( "$value" == -* && !( "$option" == "--json" && "$value" == "-" ) ) ]]; then
        echo "Error: ${option} requires a value." >&2
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idt-threshold)
            require_option_value "$1" "${2:-}"
            IDT_THRESHOLD="$2"
            shift 2
            ;;
        --idt-threshold=*)
            IDT_THRESHOLD="${1#*=}"
            require_option_value "--idt-threshold" "$IDT_THRESHOLD"
            shift
            ;;
        --json)
            require_option_value "$1" "${2:-}"
            JSON_OUTPUT="$2"
            shift 2
            ;;
        --json=*)
            JSON_OUTPUT="${1#*=}"
            require_option_value "--json" "$JSON_OUTPUT"
            shift
            ;;
        --use-filtered)
            USE_FILTERED=1
            shift
            ;;
        -o|--report-output)
            require_option_value "$1" "${2:-}"
            REPORT_OUTPUT="$2"
            shift 2
            ;;
        --report-output=*)
            REPORT_OUTPUT="${1#*=}"
            require_option_value "--report-output" "$REPORT_OUTPUT"
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
        --collect-only)
            COLLECT_ONLY=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            POSITIONAL_ARGS+=("$@")
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL_ARGS[@]+"${POSITIONAL_ARGS[@]}"}"

if [[ $# -ne 1 ]]; then
    if [[ $# -lt 1 ]]; then
        echo "Error: input_dir is required." >&2
    else
        echo "Error: exactly one input_dir is required." >&2
    fi
    usage >&2
    exit 1
fi

INPUT_DIR="$1"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input_dir not found or not a directory: $INPUT_DIR" >&2
    exit 1
fi

if [[ ! -f "$COLLECT_PY" ]]; then
    echo "Error: collect_run_summary.py not found at $COLLECT_PY" >&2
    exit 1
fi

if [[ ! -f "$CONTRACT_PY" ]]; then
    echo "Error: report_contract.py not found at $CONTRACT_PY" >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is not available in PATH." >&2
    exit 1
fi

if [[ "$JSON_OUTPUT" == "-" && "$COLLECT_ONLY" -ne 1 ]]; then
    echo "Error: --json - can only be used with --collect-only." >&2
    exit 1
fi

python3 "$CONTRACT_PY" validate-format "$OUTPUT_FORMAT"

# --- Step 1: collect_run_summary.py ---

COLLECT_ARGS=("$INPUT_DIR")
if [[ -n "$IDT_THRESHOLD" ]]; then
    COLLECT_ARGS+=("--idt-threshold" "$IDT_THRESHOLD")
fi
if [[ -n "$JSON_OUTPUT" ]]; then
    COLLECT_ARGS+=("--output" "$JSON_OUTPUT")
fi

echo "--- Step 1: collect_run_summary.py ---"
python3 "$COLLECT_PY" "${COLLECT_ARGS[@]}"

if [[ "$COLLECT_ONLY" -eq 1 ]]; then
    echo "--- collect-only mode: done ---"
    exit 0
fi

# --- Step 2: derive JSON path and run render_report.sh ---

SUMMARY_PATH_ARGS=("summary-path" "$INPUT_DIR")
if [[ -n "$JSON_OUTPUT" ]]; then
    SUMMARY_PATH_ARGS+=("--output" "$JSON_OUTPUT")
fi
if [[ "$USE_FILTERED" -eq 1 ]]; then
    SUMMARY_PATH_ARGS+=("--filtered")
fi
REPORT_JSON="$(python3 "$CONTRACT_PY" "${SUMMARY_PATH_ARGS[@]}")"

if [[ ! -f "$REPORT_JSON" ]]; then
    echo "Error: Expected JSON not found: $REPORT_JSON" >&2
    exit 1
fi

if [[ ! -f "$RENDER_SH" ]]; then
    echo "Error: render_report.sh not found at $RENDER_SH" >&2
    exit 1
fi

RENDER_ARGS=("-j" "$REPORT_JSON" "-f" "$OUTPUT_FORMAT")
if [[ -n "$REPORT_OUTPUT" ]]; then
    RENDER_ARGS+=("-o" "$REPORT_OUTPUT")
fi

echo "--- Step 2: render_report.sh (json: $(basename "$REPORT_JSON")) ---"
bash "$RENDER_SH" "${RENDER_ARGS[@]}"
