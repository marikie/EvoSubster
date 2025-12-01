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
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

IDENTITY_PREFIX = "# substitution percent identity:"


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
        help="Optional file path where the JSON summary will be written. Defaults to stdout.",
    )
    parser.add_argument(
        "--ncbi-cli",
        default="datasets",
        help="NCBI CLI executable used to fetch metadata (default: datasets).",
    )
    parser.add_argument(
        "--ncbi-timeout",
        type=int,
        default=30,
        help="Timeout in seconds for the NCBI CLI invocation.",
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


def matches_species_suffix(path: Path, species_suffix: str) -> bool:
    parts = path.stem.split("_")
    for part in parts:
        if part.endswith(species_suffix) and not part.isdigit():
            return True
    return False


def filter_for_species(paths: Iterable[Path], species_suffix: str) -> List[Path]:
    filtered = [path for path in paths if matches_species_suffix(path, species_suffix)]
    return sorted(filtered)


def select_first_for_species(run_dir: Path, pattern: str, species_suffix: str) -> Optional[Path]:
    matches = filter_for_species(run_dir.glob(pattern), species_suffix)
    if not matches:
        return None
    return matches[0]


def select_all_for_species(run_dir: Path, pattern: str, species_suffix: str) -> List[Path]:
    return filter_for_species(run_dir.glob(pattern), species_suffix)


def extract_accession_from_filename(path: Path) -> Optional[str]:
    parts = path.stem.split("_")
    if len(parts) < 2:
        return None
    return "_".join(parts[:2])


def load_json_file(path: Path) -> Tuple[Optional[dict], Optional[str]]:
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except json.JSONDecodeError as exc:
        return None, f"Failed to parse JSON at {path}: {exc}"
    except OSError as exc:
        return None, f"Failed to read {path}: {exc}"


def fetch_ncbi_metadata(
    accession: str, cache_dir: Path, cli: str, timeout: int
) -> Dict[str, Optional[str]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    json_path = cache_dir / f"{accession}.json"
    command = [cli, "summary", "genome", "accession", accession]

    result: Dict[str, Optional[str]] = {
        "species_name": None,
        "command": " ".join(command),
        "returncode": None,
        "stderr": None,
        "error": None,
        "cached": None,
        "json_path": str(json_path),
    }

    payload: Optional[dict] = None

    if json_path.exists():
        payload, error = load_json_file(json_path)
        if error is None and payload is not None:
            result["cached"] = True
        else:
            result["cached"] = False
            result["error"] = error

    if payload is None:
        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except FileNotFoundError:
            result["error"] = f"CLI not found: {cli}"
            return result
        except subprocess.TimeoutExpired:
            result["error"] = f"CLI timed out after {timeout} seconds"
            return result

        result["returncode"] = str(completed.returncode)
        stderr = completed.stderr.strip() if completed.stderr else None
        if stderr:
            result["stderr"] = stderr

        if completed.returncode != 0:
            result["error"] = "CLI exited with a non-zero status"
            return result

        stdout = completed.stdout.strip()
        if not stdout:
            result["error"] = "CLI returned an empty response"
            return result

        try:
            json_path.write_text(stdout, encoding="utf-8")
        except OSError as exc:
            result["error"] = f"Failed to write metadata JSON: {exc}"
            return result

        payload, error = load_json_file(json_path)
        if error:
            result["error"] = error
            return result
        result["cached"] = False
        result["error"] = None

    reports = payload.get("reports") if payload else None
    if not reports:
        result["error"] = "No reports returned in metadata payload"
        return result

    organism = reports[0].get("organism", {})
    organism_name = organism.get("organism_name")
    if not organism_name:
        result["error"] = "Organism name missing in metadata payload"
        return result

    result["species_name"] = organism_name
    return result


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
    return [line.strip() for line in content if line.strip()]


def process_species(
    run_dir: Path,
    species_suffix: str,
    sbst_ratio_path: Optional[Path],
    ratio_line_index: int,
    cli: str,
    timeout: int,
) -> Tuple[Dict[str, object], List[str]]:
    data: Dict[str, object] = {}
    issues: List[str] = []

    tsv_path = select_first_for_species(run_dir, f"*_*{species_suffix}_*.tsv", species_suffix)
    if tsv_path:
        data["tsv_file"] = str(tsv_path)
        accession = extract_accession_from_filename(tsv_path)
        data["accession"] = accession
        if accession:
            metadata = fetch_ncbi_metadata(accession, tsv_path.parent, cli, timeout)
            data["ncbi_lookup"] = metadata
            if metadata.get("error"):
                issues.append(
                    f"NCBI metadata lookup failed for {accession}: {metadata['error']}"
                )
        else:
            issues.append(f"Failed to derive accession from {tsv_path.name}")
    else:
        issues.append(f"No TSV matching '*_*{species_suffix}_*.tsv' found.")

    train_path = select_first_for_species(run_dir, f"*1*{species_suffix}_*.train", species_suffix)
    if train_path:
        data["train_file"] = str(train_path)
        identity_line = extract_identity_line(train_path)
        if identity_line:
            data["train_identity_line"] = identity_line
        else:
            issues.append(
                f"No '{IDENTITY_PREFIX}' line available in {train_path.name}."
            )
    else:
        issues.append(f"No train file matching '*1*{species_suffix}_*.train' found.")

    if sbst_ratio_path:
        ratio_value = extract_sbst_ratio_line(sbst_ratio_path, ratio_line_index)
        if ratio_value:
            data["sbst_ratio_entry"] = ratio_value
        else:
            issues.append(
                f"Could not extract substitution ratio line {ratio_line_index} from "
                f"{sbst_ratio_path.name}."
            )
    else:
        issues.append("Missing sbstRatio*_maflinked.out file.")

    gc_path = select_first_for_species(run_dir, f"*{species_suffix}_gcContent*.out", species_suffix)
    if gc_path:
        data["gc_content_file"] = str(gc_path)
        data["gc_content"] = read_file_lines(gc_path)
    else:
        issues.append(f"No GC content file matching '*{species_suffix}_gcContent*.out' found.")

    data["pdfs"] = {
        "maflinked_norm": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_norm.pdf", species_suffix
            )
        ],
        "maflinked_logratio": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_logRatio*.pdf", species_suffix
            )
        ],
        "maflinked_ncds_norm": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_ncds_norm.pdf", species_suffix
            )
        ],
        "maflinked_ncds_logratio": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_ncds_logRatio*.pdf", species_suffix
            )
        ],
        "maflinked_dinuc_tsv": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_dinuc.tsv.pdf", species_suffix
            )
        ],
        "maflinked_dinuc_ncds_tsv": [
            str(path)
            for path in select_all_for_species(
                run_dir, f"*_*{species_suffix}_*_maflinked_dinuc_ncds.tsv.pdf", species_suffix
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
                f"Missing required PDF matching '*_*{species_suffix}_{description}'."
            )

    return data, issues


def process_dataset(
    dataset_dir: Path, cli: str, timeout: int
) -> Tuple[Optional[Dict[str, object]], List[str]]:
    dataset_data: Dict[str, object] = {"dataset": dataset_dir.name}
    issues: List[str] = []

    latest_run = find_latest_run_dir(dataset_dir)
    if not latest_run:
        issues.append(f"No dated run directories found under {dataset_dir}.")
        return None, issues

    dataset_data["latest_run"] = latest_run.name
    dataset_data["run_dir"] = str(latest_run)

    sbst_ratio_path = select_first(latest_run, "sbstRatio*_maflinked.out")

    species2_data, species2_issues = process_species(
        latest_run, "2", sbst_ratio_path, 1, cli, timeout
    )
    species3_data, species3_issues = process_species(
        latest_run, "3", sbst_ratio_path, 2, cli, timeout
    )

    dataset_data["species2"] = species2_data
    dataset_data["species3"] = species3_data
    issues.extend(species2_issues)
    issues.extend(species3_issues)

    return dataset_data, issues


def build_summary(root: Path, cli: str, timeout: int) -> Dict[str, object]:
    summary: Dict[str, object] = {"input_root": str(root), "datasets": [], "issues": []}
    dataset_entries = []
    aggregated_issues: List[str] = []

    processed_candidate = False

    for entry in sorted(root.iterdir()):
        if not entry.is_dir() or entry.name.isdigit():
            continue
        processed_candidate = True
        dataset_data, issues = process_dataset(entry, cli, timeout)
        if dataset_data:
            dataset_entries.append(dataset_data)
        if issues:
            aggregated_issues.extend([f"{entry.name}: {msg}" for msg in issues])

    if not processed_candidate:
        dataset_data, issues = process_dataset(root, cli, timeout)
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

    summary = build_summary(input_root, args.ncbi_cli, args.ncbi_timeout)
    output_text = json.dumps(summary, ensure_ascii=False, indent=2)

    if args.output:
        output_path = Path(args.output).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(output_text, encoding="utf-8")

    print(output_text)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

