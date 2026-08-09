#!/usr/bin/env python3

"""
Collect run-level metadata and artifact paths for downstream reporting.

The script walks dataset directories under the provided root, identifies the most
recent run (by lexicographically largest directory name), gathers summary values,
retrieves organism names from the NCBI CLI, and emits a JSON payload that an
R Markdown report can consume.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

# The thesis (section 2.3) trio selection rule lives in src/select/trio_filter.py so
# that trio_selection.R and this reporter apply the identical rule.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "select"))
from trio_filter import evaluate_filter_status  # noqa: E402

IDENTITY_PREFIX = "# substitution percent identity:"
GC_CONTENT_PREFIX = "Total GC content: "


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect summary values and artifact paths from the latest run under each dataset "
            "directory, enriching the results with organism metadata from the NCBI CLI."
        )
    )
    parser.add_argument(
        "input_dir",
        help="Root directory that contains dataset subdirectories (e.g., /home/.../data/fungi).",
    )
    parser.add_argument(
        "--output",
        help=(
            "Optional file path for the JSON summary. "
            "Defaults to <input_dir>/<input_dir_name>_summary.json."
        ),
    )
    parser.add_argument(
        "--idt-threshold",
        type=float,
        default=80.0,
        help="Identity threshold; if any of idt12/idt13/idt23 is below this, the dataset is filtered out (default: 80).",
    )
    return parser.parse_args()


def find_latest_run_dir(dataset_dir: Path) -> Optional[Path]:
    run_dirs = [
        entry for entry in dataset_dir.iterdir() if entry.is_dir() and entry.name.isdigit()
    ]
    if not run_dirs:
        return None
    return max(run_dirs, key=lambda path: path.name)


def select_first(run_dir: Path, pattern: str) -> Optional[Path]:
    matches = sorted(run_dir.glob(pattern))
    if not matches:
        return None
    return matches[0]


def select_all(run_dir: Path, pattern: str) -> List[Path]:
    return sorted(run_dir.glob(pattern))


def matches_short_name(path: Path, short_name: str) -> bool:
    if not short_name:
        return False
    parts = path.stem.split("_")
    return any(part == short_name for part in parts)


def filter_for_short_name(paths: Iterable[Path], short_name: str) -> List[Path]:
    filtered = [path for path in paths if matches_short_name(path, short_name)]
    return sorted(filtered)


def select_first_for_short_name(run_dir: Path, pattern: str, short_name: str) -> Optional[Path]:
    matches = filter_for_short_name(run_dir.glob(pattern), short_name)
    if not matches:
        return None
    return matches[0]


def select_all_for_short_name(run_dir: Path, pattern: str, short_name: str) -> List[Path]:
    return filter_for_short_name(run_dir.glob(pattern), short_name)


def load_json_file(path: Path) -> Tuple[Optional[dict], Optional[str]]:
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except json.JSONDecodeError as exc:
        return None, f"Failed to parse JSON at {path}: {exc}"
    except OSError as exc:
        return None, f"Failed to read {path}: {exc}"


def read_taxonomy_lineage(metadata_json_path: str) -> Tuple[Optional[dict], Optional[str]]:
    """Read pre-fetched taxonomy classification from taxonomy_*.json alongside the metadata JSON."""
    meta_path = Path(metadata_json_path)
    tax_path = meta_path.parent / f"taxonomy_{meta_path.name}"
    if not tax_path.exists():
        return None, f"Taxonomy JSON not found: {tax_path}"
    data, err = load_json_file(tax_path)
    if err or not data:
        return None, err
    reports = data.get("reports") or []
    if not reports:
        return None, "No reports in taxonomy JSON"
    classification = (reports[0].get("taxonomy") or {}).get("classification") or {}
    ranks = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    result = {}
    for rank in ranks:
        entry = classification.get(rank)
        if entry and entry.get("name"):
            result[rank] = entry["name"]
    if not result:
        return None, "No classification data found in taxonomy JSON"
    return result, None


def load_manifest(run_dir: Path) -> Tuple[Optional[dict], Optional[str]]:
    manifest_path = run_dir / "metadata" / "metadata_manifest.json"
    if not manifest_path.exists():
        return None, f"Metadata manifest not found: {manifest_path}"
    payload, error = load_json_file(manifest_path)
    if error:
        return None, error
    payload["manifest_path"] = str(manifest_path)
    return payload, None


def build_manifest_slot_map(manifest: dict) -> Dict[str, dict]:
    slot_map: Dict[str, dict] = {}
    for entry in manifest.get("organisms", []):
        slot = entry.get("slot")
        if slot:
            slot_map[slot] = entry
    return slot_map


def sanitize_metadata(metadata: dict) -> dict:
    """Drop duplicated keys that are promoted to the species root."""
    filtered = {
        key: value
        for key, value in metadata.items()
        if key not in {"accession", "metadata_json"}
    }
    return filtered


def extract_identity_line(train_path: Path) -> Optional[str]:
    lines = [
        line.strip()
        for line in train_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(IDENTITY_PREFIX)
    ]
    if len(lines) >= 2:
        return lines[-2]
    if lines:
        return lines[0]
    return None


def extract_identity_value(identity_line: str) -> str:
    if ":" in identity_line:
        _, remainder = identity_line.split(":", 1)
        return remainder.strip()
    return identity_line.strip()


def extract_sbst_ratio_line(sbst_path: Path, line_index: int) -> Optional[str]:
    lines = sbst_path.read_text(encoding="utf-8").splitlines()
    if line_index < 1 or line_index > len(lines):
        return None
    target = lines[line_index - 1]
    if ":" in target:
        return target.split(":", 1)[1].strip()
    return target.strip()


def read_file_lines(path: Path) -> List[str]:
    content = path.read_text(encoding="utf-8").splitlines()
    normalized: List[str] = []
    for raw_line in content:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(GC_CONTENT_PREFIX):
            line = line[len(GC_CONTENT_PREFIX) :].lstrip()
        normalized.append(line)
    return normalized


def process_species(
    run_dir: Path,
    metadata: dict,
    sbst_ratio_path: Optional[Path],
    ratio_line_index: Optional[int],
    *,
    include_tsv: bool = True,
    include_sbst_ratio: bool = True,
    include_pdfs: bool = True,
) -> Tuple[Dict[str, object], List[str]]:
    data: Dict[str, object] = {"metadata": sanitize_metadata(metadata)}
    issues: List[str] = []

    slot = metadata.get("slot", "unknown")
    short_name = metadata.get("short_name")
    if not short_name:
        issues.append(f"{slot}: Metadata missing short_name.")
        return data, issues

    accession = metadata.get("accession")
    if accession:
        data["accession"] = accession
    metadata_json = metadata.get("metadata_json")
    if metadata_json:
        data["metadata_json"] = metadata_json
        lineage, tax_err = read_taxonomy_lineage(metadata_json)
        if lineage:
            data["taxonomy"] = lineage
        elif tax_err:
            issues.append(f"{slot}: {tax_err}")

    if include_tsv:
        tsv_pattern = f"statistics/{short_name}/singlenuc/*_*{short_name}_*.tsv"
        tsv_path = select_first_for_short_name(run_dir, tsv_pattern, short_name)
        if tsv_path:
            data["tsv_file"] = str(tsv_path)
        else:
            issues.append(f"{slot}: No TSV matching '{tsv_pattern}' found.")

    if include_sbst_ratio:
        if sbst_ratio_path and ratio_line_index is not None:
            ratio_value = extract_sbst_ratio_line(sbst_ratio_path, ratio_line_index)
            if ratio_value:
                data["sbst_ratio_entry"] = ratio_value
            else:
                issues.append(
                    f"{slot}: Could not extract substitution ratio line {ratio_line_index} from "
                    f"{sbst_ratio_path.name}."
                )
        else:
            issues.append(f"{slot}: Missing statistics/misc/sbstRatio*.out file.")

    gc_pattern = f"statistics/misc/*{short_name}_gcContent*.out"
    gc_path = select_first_for_short_name(run_dir, gc_pattern, short_name)
    if gc_path:
        data["gc_content_file"] = str(gc_path)
        data["gc_content"] = read_file_lines(gc_path)
    else:
        issues.append(
            f"{slot}: No GC content file matching '{gc_pattern}' found."
        )

    if include_pdfs:
        fig_sn = f"figs/{short_name}/singlenuc"
        fig_dn = f"figs/{short_name}/dinuc"
        data["pdfs"] = {
            "norm": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_sn}/ratio/*_*{short_name}_*_norm.pdf", short_name
                )
                if "_ncds_norm" not in path.name
            ],
            "logratio": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_sn}/log-ratio/*_*{short_name}_*_logRatio*.pdf", short_name
                )
                if "_ncds_logRatio" not in path.name
            ],
            "ncds_norm": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_sn}/ratio/*_*{short_name}_*_ncds_norm.pdf", short_name
                )
            ],
            "ncds_logratio": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_sn}/log-ratio/*_*{short_name}_*_ncds_logRatio*.pdf", short_name
                )
            ],
            "dinuc_tsv": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_dn}/*_*{short_name}_*_dinuc.tsv.pdf", short_name
                )
            ],
            "dinuc_ncds_tsv": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"{fig_dn}/*_*{short_name}_*_dinuc_ncds.tsv.pdf", short_name
                )
            ],
        }

        required_pdfs = {
            "norm": "norm.pdf",
            "logratio": "logRatio*.pdf",
            "dinuc_tsv": "dinuc.tsv.pdf",
        }
        for key, description in required_pdfs.items():
            if not data["pdfs"][key]:
                issues.append(
                    f"{slot}: Missing required PDF matching '*_*{short_name}_{description}'."
                )

    return data, issues


def collect_identity_metrics(
    run_dir: Path, slot_map: Dict[str, dict]
) -> Tuple[Dict[str, str], List[str]]:
    identity: Dict[str, str] = {}
    issues: List[str] = []

    def short_name_for(slot: str) -> Optional[str]:
        entry = slot_map.get(slot)
        return entry.get("short_name") if entry else None

    pair_specs = [
        ("idt_12", "org1", "org2"),
        ("idt_13", "org1", "org3"),
        ("idt_23", "org2", "org3"),
    ]

    for key, slot_a, slot_b in pair_specs:
        short_a = short_name_for(slot_a)
        short_b = short_name_for(slot_b)
        if not short_a or not short_b:
            issues.append(f"{key}: Missing short_name for {slot_a} or {slot_b}.")
            continue

        pattern = f"intermediateFiles/*{short_a}2{short_b}_*.train"
        train_path = select_first(run_dir, pattern)
        if not train_path:
            issues.append(f"{key}: No train file matching '{pattern}' found.")
            continue

        identity_line = extract_identity_line(train_path)
        if not identity_line:
            issues.append(f"{key}: '{IDENTITY_PREFIX}' line missing in {train_path.name}.")
            continue

        raw_value = extract_identity_value(identity_line)
        value_with_unit = raw_value if raw_value.endswith("%") else f"{raw_value} %"
        identity[key] = value_with_unit

    return identity, issues


def process_dataset(dataset_dir: Path, idt_threshold: float) -> Tuple[Optional[Dict[str, object]], List[str]]:
    dataset_data: Dict[str, object] = {"dataset": dataset_dir.name}
    issues: List[str] = []

    latest_run = find_latest_run_dir(dataset_dir)
    if not latest_run:
        issues.append(f"No dated run directories found under {dataset_dir}.")
        return None, issues

    dataset_data["latest_run"] = latest_run.name
    dataset_data["run_dir"] = str(latest_run)

    manifest, manifest_error = load_manifest(latest_run)
    if manifest_error or not manifest:
        issues.append(manifest_error or "Unknown metadata manifest error.")
        return None, issues
    dataset_data["metadata_manifest"] = manifest.get("manifest_path")
    slot_map = build_manifest_slot_map(manifest)

    sbst_ratio_path = select_first(latest_run, "statistics/misc/sbstRatio*.out")

    species1_meta = slot_map.get("org1")
    if not species1_meta:
        species1_data = {}
        issues.append("org1: Missing metadata entry in manifest.")
    else:
        species1_data, species1_issues = process_species(
            latest_run,
            species1_meta,
            None,
            None,
            include_tsv=False,
            include_sbst_ratio=False,
            include_pdfs=False,
        )
        issues.extend(species1_issues)
    dataset_data["species1"] = species1_data

    species2_meta = slot_map.get("org2")
    if not species2_meta:
        species2_data = {}
        issues.append("org2: Missing metadata entry in manifest.")
    else:
        species2_data, species2_issues = process_species(
            latest_run, species2_meta, sbst_ratio_path, 1
        )
        issues.extend(species2_issues)
    dataset_data["species2"] = species2_data

    species3_meta = slot_map.get("org3")
    if not species3_meta:
        species3_data = {}
        issues.append("org3: Missing metadata entry in manifest.")
    else:
        species3_data, species3_issues = process_species(
            latest_run, species3_meta, sbst_ratio_path, 2
        )
        issues.extend(species3_issues)
    dataset_data["species3"] = species3_data

    identity_values, identity_issues = collect_identity_metrics(latest_run, slot_map)
    dataset_data.update(identity_values)
    issues.extend(identity_issues)

    filter_status, filter_issues = evaluate_filter_status(slot_map, identity_values, idt_threshold)
    dataset_data["filter_status"] = filter_status
    issues.extend(filter_issues)

    return dataset_data, issues


def format_issue_list(issues_map: Dict[str, List[str]], allowed: Optional[Set[str]] = None) -> List[str]:
    messages: List[str] = []
    for dataset_name in sorted(issues_map.keys()):
        if allowed is not None and dataset_name not in allowed:
            continue
        for issue in issues_map[dataset_name]:
            messages.append(f"{dataset_name}: {issue}")
    return messages


def build_summary(root: Path, idt_threshold: float) -> Tuple[Dict[str, object], Dict[str, object]]:
    summary_all: Dict[str, object] = {"input_root": str(root), "datasets": [], "issues": []}
    summary_filtered: Dict[str, object] = {"input_root": str(root), "datasets": [], "issues": []}
    dataset_entries: List[Dict[str, object]] = []
    filtered_entries: List[Dict[str, object]] = []
    aggregated_issues: Dict[str, List[str]] = {}

    processed_candidate = False

    for entry in sorted(root.iterdir()):
        if not entry.is_dir() or entry.name.isdigit():
            continue
        processed_candidate = True
        dataset_data, issues = process_dataset(entry, idt_threshold)
        dataset_name = entry.name
        if dataset_data:
            dataset_entries.append(dataset_data)
            filter_info = dataset_data.get("filter_status", {})
            if not filter_info.get("excluded"):
                filtered_entries.append(dataset_data)
            dataset_name = dataset_data.get("dataset", dataset_name)
        if issues:
            aggregated_issues.setdefault(dataset_name, []).extend(issues)

    if not processed_candidate:
        dataset_data, issues = process_dataset(root, idt_threshold)
        dataset_name = root.name
        if dataset_data:
            dataset_entries.append(dataset_data)
            filter_info = dataset_data.get("filter_status", {})
            if not filter_info.get("excluded"):
                filtered_entries.append(dataset_data)
            dataset_name = dataset_data.get("dataset", dataset_name)
        if issues:
            aggregated_issues.setdefault(dataset_name, []).extend(issues)

    summary_all["datasets"] = dataset_entries
    summary_all["issues"] = format_issue_list(aggregated_issues)

    filtered_names = {entry.get("dataset") for entry in filtered_entries if entry.get("dataset")}
    summary_filtered["datasets"] = filtered_entries
    summary_filtered["issues"] = format_issue_list(aggregated_issues, filtered_names)
    return summary_all, summary_filtered


def main() -> None:
    args = parse_args()
    input_root = Path(args.input_dir).resolve()

    if not input_root.is_dir():
        print(f"Error: input directory does not exist: {input_root}", file=sys.stderr)
        sys.exit(1)

    summary_all, summary_filtered = build_summary(input_root, args.idt_threshold)
    output_text = json.dumps(summary_all, ensure_ascii=False, indent=2)
    filtered_output_text = json.dumps(summary_filtered, ensure_ascii=False, indent=2)

    output_path: Optional[Path]
    filtered_output_path: Optional[Path]

    def derive_filtered_path(base_path: Path) -> Path:
        suffix = base_path.suffix or ".json"
        stem = base_path.stem if base_path.suffix else base_path.name
        return base_path.with_name(f"{stem}_filtered{suffix}")

    if args.output:
        if args.output.strip() == "-":
            output_path = None
            filtered_output_path = None
        else:
            output_path = Path(args.output).resolve()
            filtered_output_path = derive_filtered_path(output_path)
    else:
        default_name = f"{input_root.name}_summary.json"
        output_path = input_root / default_name
        filtered_output_path = input_root / f"{input_root.name}_summary_filtered.json"

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(output_text, encoding="utf-8")

    if filtered_output_path is not None:
        filtered_output_path.parent.mkdir(parents=True, exist_ok=True)
        filtered_output_path.write_text(filtered_output_text, encoding="utf-8")

    print(output_text)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

