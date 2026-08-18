#!/usr/bin/env python3

"""Shared path and JSON-contract utilities for the report pipeline."""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REQUIRED_DATASET_KEYS = (
    "species1",
    "species2",
    "species3",
    "idt_12",
    "idt_13",
    "idt_23",
)

FORMAT_EXTENSIONS = {
    "word_document": ".docx",
    "pdf_document": ".pdf",
    "html_document": ".html",
}


def validate_output_format(output_format: str) -> None:
    if output_format not in FORMAT_EXTENSIONS:
        supported = ", ".join(FORMAT_EXTENSIONS)
        raise ValueError(
            f"unsupported output format '{output_format}'; expected one of: {supported}"
        )


def derive_filtered_path(summary_path: Path) -> Path:
    suffix = summary_path.suffix or ".json"
    stem = summary_path.stem if summary_path.suffix else summary_path.name
    return summary_path.with_name(f"{stem}_filtered{suffix}")


def summary_paths(input_root: Path, output: Optional[str] = None) -> Tuple[Path, Path]:
    root = input_root.resolve()
    summary_path = (
        Path(output).resolve()
        if output
        else root / f"{root.name}_summary.json"
    )
    return summary_path, derive_filtered_path(summary_path)


def report_output_path(
    summary_path: Path, output_format: str, output: Optional[str] = None
) -> Path:
    validate_output_format(output_format)
    if output:
        return Path(output).resolve()
    extension = FORMAT_EXTENSIONS[output_format]
    return summary_path.resolve().with_suffix(extension)


def preview_dir(summary_path: Path) -> Path:
    path = summary_path.resolve().parent / "tmp"
    path.mkdir(parents=True, exist_ok=True)
    return path


def validate_summary(payload: object) -> Tuple[List[str], List[str]]:
    if not isinstance(payload, dict):
        return [], ["summary JSON must contain an object at the top level."]

    datasets = payload.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        return [], ["summary JSON contains no datasets."]

    names: List[str] = []
    errors: List[str] = []
    for entry in datasets:
        if not isinstance(entry, dict):
            errors.append("summary JSON contains a dataset entry that is not an object.")
            continue
        name = entry.get("dataset") or "unknown dataset"
        names.append(str(name))
        missing = [key for key in REQUIRED_DATASET_KEYS if key not in entry]
        if missing:
            errors.append(
                f"{name} is missing required fields from collect_run_summary.py: "
                f"{', '.join(missing)}"
            )
            continue

        for species_key in ("species1", "species2", "species3"):
            species = entry[species_key]
            if not isinstance(species, dict):
                errors.append(f"{name}.{species_key} must be an object.")
                continue
            for object_key in ("metadata", "taxonomy"):
                value = species.get(object_key)
                if value is not None and not isinstance(value, dict):
                    errors.append(
                        f"{name}.{species_key}.{object_key} must be an object."
                    )
            gc_content = species.get("gc_content")
            if gc_content is not None and (
                not isinstance(gc_content, list)
                or not all(isinstance(value, str) for value in gc_content)
            ):
                errors.append(
                    f"{name}.{species_key}.gc_content must be an array of strings."
                )
            pdfs = species.get("pdfs")
            if pdfs is not None:
                if not isinstance(pdfs, dict):
                    errors.append(f"{name}.{species_key}.pdfs must be an object.")
                else:
                    for pdf_key, paths in pdfs.items():
                        if not isinstance(paths, list) or not all(
                            isinstance(path, str) for path in paths
                        ):
                            errors.append(
                                f"{name}.{species_key}.pdfs.{pdf_key} must be an array of strings."
                            )

        for identity_key in ("idt_12", "idt_13", "idt_23"):
            if not isinstance(entry[identity_key], str):
                errors.append(f"{name}.{identity_key} must be a string.")
    return names, errors


def load_and_validate_summary(path: Path) -> List[str]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"failed to parse summary JSON at {path}: {exc}") from exc
    except OSError as exc:
        raise ValueError(f"failed to read summary JSON at {path}: {exc}") from exc

    names, errors = validate_summary(payload)
    if errors:
        raise ValueError("\n".join(errors))
    return names


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    summary_parser = subparsers.add_parser("summary-path")
    summary_parser.add_argument("input_dir")
    summary_parser.add_argument("--output")
    summary_parser.add_argument("--filtered", action="store_true")

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("json_path")

    report_parser = subparsers.add_parser("report-path")
    report_parser.add_argument("json_path")
    report_parser.add_argument("output_format")
    report_parser.add_argument("--output")

    preview_parser = subparsers.add_parser("preview-dir")
    preview_parser.add_argument("json_path")

    format_parser = subparsers.add_parser("validate-format")
    format_parser.add_argument("output_format")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "summary-path":
            summary_path, filtered_path = summary_paths(
                Path(args.input_dir), args.output
            )
            print(filtered_path if args.filtered else summary_path)
        elif args.command == "validate":
            print(",".join(load_and_validate_summary(Path(args.json_path))))
        elif args.command == "report-path":
            print(
                report_output_path(
                    Path(args.json_path), args.output_format, args.output
                )
            )
        elif args.command == "preview-dir":
            print(preview_dir(Path(args.json_path)))
        elif args.command == "validate-format":
            validate_output_format(args.output_format)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
