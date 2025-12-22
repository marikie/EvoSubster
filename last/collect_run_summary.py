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
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

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


def extract_genus_label(metadata_entry: dict) -> Optional[str]:
    candidates = [
        metadata_entry.get("raw_organism_name"),
        metadata_entry.get("ncbi_full_name"),
        metadata_entry.get("directory_name"),
        metadata_entry.get("short_name"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        sanitized = candidate.replace("_", " ").strip()
        if not sanitized:
            continue
        token = sanitized.split()[0]
        if token:
            return token
    return None


def evaluate_genus_condition(slot_map: Dict[str, dict]) -> Tuple[Optional[bool], List[str]]:
    patterns, issues = evaluate_genus_patterns(slot_map)
    if patterns is None:
        return None, issues
    return patterns["pattern1"], issues


def evaluate_genus_patterns(
    slot_map: Dict[str, dict]
) -> Tuple[Optional[Dict[str, bool]], List[str]]:
    issues: List[str] = []
    required_slots = ("org1", "org2", "org3")
    genus_values: Dict[str, str] = {}

    for slot in required_slots:
        entry = slot_map.get(slot)
        if not entry:
            issues.append(f"{slot}: Missing metadata entry for genus evaluation.")
            return None, issues
        genus_label = extract_genus_label(entry)
        if not genus_label:
            issues.append(f"{slot}: Unable to derive genus from metadata.")
            return None, issues
        genus_values[slot] = genus_label.lower()

    pattern1 = (
        genus_values["org1"] != genus_values["org2"]
        and genus_values["org1"] != genus_values["org3"]
        and genus_values["org2"] == genus_values["org3"]
    )
    pattern3 = (
        genus_values["org1"] == genus_values["org2"]
        and genus_values["org2"] == genus_values["org3"]
    )
    pattern4 = (
        genus_values["org1"] != genus_values["org2"]
        and genus_values["org1"] != genus_values["org3"]
        and genus_values["org2"] != genus_values["org3"]
    )
    pattern5 = (
        (
            genus_values["org1"] != genus_values["org2"]
            and genus_values["org1"] == genus_values["org3"]
            and genus_values["org2"] != genus_values["org3"]
        )
        or (
            genus_values["org1"] == genus_values["org2"]
            and genus_values["org1"] != genus_values["org3"]
            and genus_values["org2"] != genus_values["org3"]
        )
    )

    patterns = {
        "pattern1": pattern1,
        "pattern3": pattern3,
        "pattern4": pattern4,
        "pattern5": pattern5,
    }
    return patterns, issues


def parse_identity_percentage(value: str) -> Optional[float]:
    if not value:
        return None
    sanitized = value.replace("%", " ")
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", sanitized)
    if not match:
        return None
    try:
        return float(match.group())
    except ValueError:
        return None


def evaluate_identity_condition(
    identity_values: Dict[str, str]
) -> Tuple[Optional[bool], List[str], Optional[Dict[str, float]]]:
    issues: List[str] = []
    numeric_values: Dict[str, float] = {}
    required_keys = ("idt_12", "idt_13", "idt_23")

    for key in required_keys:
        raw_value = identity_values.get(key)
        if raw_value is None:
            issues.append(f"{key}: Identity value missing; cannot evaluate condition.")
            return None, issues, None
        parsed = parse_identity_percentage(raw_value)
        if parsed is None:
            issues.append(f"{key}: Unable to parse numeric identity from '{raw_value}'.")
            return None, issues, None
        numeric_values[key] = parsed

    condition_met = numeric_values["idt_12"] < numeric_values["idt_23"] and numeric_values["idt_13"] < numeric_values["idt_23"]
    return condition_met, issues, numeric_values


def evaluate_filter_status(
    slot_map: Dict[str, dict],
    identity_values: Dict[str, str],
    idt_threshold: float,
) -> Tuple[Dict[str, Optional[bool]], List[str]]:
    issues: List[str] = []
    genus_patterns, genus_issues = evaluate_genus_patterns(slot_map)
    identity_condition, identity_issues, numeric_identity_values = evaluate_identity_condition(
        identity_values
    )
    issues.extend(genus_issues)
    issues.extend(identity_issues)

    genus_condition = genus_patterns["pattern1"] if genus_patterns else None
    genus_pattern_same_all = genus_patterns["pattern3"] if genus_patterns else None
    genus_pattern_all_different = genus_patterns["pattern4"] if genus_patterns else None
    genus_pattern_two_vs_one = genus_patterns["pattern5"] if genus_patterns else None

    threshold_condition: Optional[bool] = None
    if numeric_identity_values is not None:
        threshold_condition = all(value > idt_threshold for value in numeric_identity_values.values())

    filter_flag: Optional[bool] = None

    # Threshold gate: if below threshold, filter out before other rules
    if threshold_condition is False:
        filter_flag = True
    else:
        if genus_condition is True:
            filter_flag = False
        elif identity_condition is True:
            if genus_pattern_same_all is True:
                filter_flag = False
            elif genus_pattern_all_different is True:
                filter_flag = False
            elif genus_pattern_two_vs_one is True:
                filter_flag = True
        elif identity_condition is False:
            filter_flag = True

    excluded = filter_flag is True
    filter_details = {
        "genus_condition_met": genus_condition,
        "genus_pattern_same_all": genus_pattern_same_all,
        "genus_pattern_all_different": genus_pattern_all_different,
        "genus_pattern_two_vs_one": genus_pattern_two_vs_one,
        "identity_condition_met": identity_condition,
        "idt_threshold_condition": threshold_condition,
        "idt_threshold_value": idt_threshold,
        "filter_flag": filter_flag,
        "excluded": excluded,
    }
    return filter_details, issues


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

    if include_tsv:
        tsv_path = select_first_for_short_name(run_dir, f"*_*{short_name}_*.tsv", short_name)
        if tsv_path:
            data["tsv_file"] = str(tsv_path)
        else:
            issues.append(f"{slot}: No TSV matching '*_*{short_name}_*.tsv' found.")

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
            issues.append(f"{slot}: Missing sbstRatio*_maflinked.out file.")

    gc_path = select_first_for_short_name(run_dir, f"*{short_name}_gcContent*.out", short_name)
    if gc_path:
        data["gc_content_file"] = str(gc_path)
        data["gc_content"] = read_file_lines(gc_path)
    else:
        issues.append(
            f"{slot}: No GC content file matching '*{short_name}_gcContent*.out' found."
        )

    if include_pdfs:
        data["pdfs"] = {
            "maflinked_norm": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_norm.pdf", short_name
                )
            ],
            "maflinked_logratio": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_logRatio*.pdf", short_name
                )
            ],
            "maflinked_ncds_norm": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_ncds_norm.pdf", short_name
                )
            ],
            "maflinked_ncds_logratio": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_ncds_logRatio*.pdf", short_name
                )
            ],
            "maflinked_dinuc_tsv": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_dinuc.tsv.pdf", short_name
                )
            ],
            "maflinked_dinuc_ncds_tsv": [
                str(path)
                for path in select_all_for_short_name(
                    run_dir, f"*_*{short_name}_*_maflinked_dinuc_ncds.tsv.pdf", short_name
                )
            ],
        }

        required_pdfs = {
            "maflinked_norm": "maflinked_norm.pdf",
            "maflinked_logratio": "maflinked_logRatio*.pdf",
            "maflinked_dinuc_tsv": "maflinked_dinuc.tsv.pdf",
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

        pattern = f"*{short_a}2{short_b}_*.train"
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

    sbst_ratio_path = select_first(latest_run, "sbstRatio*_maflinked.out")

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

