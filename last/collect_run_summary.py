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
from typing import Dict, Iterable, List, Optional, Tuple

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


def process_dataset(dataset_dir: Path) -> Tuple[Optional[Dict[str, object]], List[str]]:
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

    return dataset_data, issues


def build_summary(root: Path) -> Dict[str, object]:
    summary: Dict[str, object] = {"input_root": str(root), "datasets": [], "issues": []}
    dataset_entries = []
    aggregated_issues: List[str] = []

    processed_candidate = False

    for entry in sorted(root.iterdir()):
        if not entry.is_dir() or entry.name.isdigit():
            continue
        processed_candidate = True
        dataset_data, issues = process_dataset(entry)
        if dataset_data:
            dataset_entries.append(dataset_data)
        if issues:
            aggregated_issues.extend([f"{entry.name}: {msg}" for msg in issues])

    if not processed_candidate:
        dataset_data, issues = process_dataset(root)
        if dataset_data:
            dataset_entries.append(dataset_data)
        if issues:
            aggregated_issues.extend([f"{root.name}: {msg}" for msg in issues])

    summary["datasets"] = dataset_entries
    summary["issues"] = aggregated_issues
    return summary


def main() -> None:
    args = parse_args()
    input_root = Path(args.input_dir).resolve()

    if not input_root.is_dir():
        print(f"Error: input directory does not exist: {input_root}", file=sys.stderr)
        sys.exit(1)

    summary = build_summary(input_root)
    output_text = json.dumps(summary, ensure_ascii=False, indent=2)

    output_path: Optional[Path]
    if args.output:
        if args.output.strip() == "-":
            output_path = None
        else:
            output_path = Path(args.output).resolve()
    else:
        default_name = f"{input_root.name}_summary.json"
        output_path = input_root / default_name

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(output_text, encoding="utf-8")

    print(output_text)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

