#!/usr/bin/env python3

"""
Summarize all trios under a lineage directory into a Markdown table.

Usage:
    python src/report/summarize_trios.py <input_dir> [--output PATH]

Each row in the output table corresponds to one organism (three rows per trio).
Columns: Trio, Non-coding, Organism name, Accession, Role, Phylum, Class,
         Order, Family, Genus, Species, Notes.

If a taxonomy JSON is missing, the script attempts to download it via
'datasets summary taxonomy taxon <tax_id>' before falling back to '—'.
"""

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


TAXONOMY_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
MISSING = "—"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

@dataclass
class OrgRow:
    trio_id: str
    has_noncoding: bool
    organism_name: str
    accession: str
    role: str           # "outgroup" or "ingroup"
    kingdom: str = MISSING
    phylum: str = MISSING
    klass: str = MISSING   # "class" is a Python keyword
    order: str = MISSING
    family: str = MISSING
    genus: str = MISSING
    species: str = MISSING
    notes: str = ""


# ---------------------------------------------------------------------------
# Directory helpers
# ---------------------------------------------------------------------------

def find_trio_dirs(root: Path) -> List[Path]:
    """Return subdirectories that look like trio combos (name has 2+ underscores)."""
    dirs = []
    for entry in sorted(root.iterdir()):
        if not entry.is_dir():
            continue
        parts = entry.name.split("_")
        if len(parts) >= 3:
            dirs.append(entry)
    return dirs


def find_latest_run_dir(trio_dir: Path) -> Optional[Path]:
    """Return the subdirectory with the lexicographically largest all-digit name."""
    run_dirs = [e for e in trio_dir.iterdir() if e.is_dir() and e.name.isdigit()]
    if not run_dirs:
        return None
    return max(run_dirs, key=lambda p: p.name)


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

def load_json_file(path: Path) -> Tuple[Optional[dict], Optional[str]]:
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except json.JSONDecodeError as exc:
        return None, f"JSON parse error at {path}: {exc}"
    except OSError as exc:
        return None, f"Cannot read {path}: {exc}"


def load_manifest(run_dir: Path) -> Tuple[Optional[dict], Optional[str]]:
    manifest_path = run_dir / "metadata" / "metadata_manifest.json"
    if not manifest_path.exists():
        return None, f"Manifest not found: {manifest_path}"
    return load_json_file(manifest_path)


# ---------------------------------------------------------------------------
# Non-coding detection
# ---------------------------------------------------------------------------

def has_noncoding_data(run_dir: Path) -> bool:
    """Return True if any *_ncds.tsv file exists directly under run_dir."""
    return bool(list(run_dir.glob("*_ncds.tsv")))


# ---------------------------------------------------------------------------
# Taxonomy helpers
# ---------------------------------------------------------------------------

def get_tax_id_from_metadata(short_name: str, accession: str, metadata_dir: Path) -> Optional[str]:
    """Extract tax_id from the per-organism NCBI metadata JSON."""
    meta_path = metadata_dir / f"{short_name}_{accession}.json"
    data, err = load_json_file(meta_path)
    if err or not data:
        return None
    reports = data.get("reports") or []
    if not reports:
        return None
    organism = reports[0].get("organism") or {}
    tax_id = organism.get("tax_id")
    if tax_id is not None:
        return str(tax_id)
    # Fallback: biosample path
    biosample = reports[0].get("biosample") or {}
    desc_organism = (biosample.get("description") or {}).get("organism") or {}
    tax_id = desc_organism.get("tax_id")
    return str(tax_id) if tax_id is not None else None


def download_taxonomy_json(tax_id: str, output_path: Path) -> Optional[str]:
    """Download taxonomy via NCBI Datasets CLI and write to output_path."""
    result = subprocess.run(
        ["datasets", "summary", "taxonomy", "taxon", tax_id],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return f"datasets taxonomy download failed for tax_id={tax_id}: {result.stderr.strip()}"
    try:
        output_path.write_text(result.stdout, encoding="utf-8")
    except OSError as exc:
        return f"Cannot write taxonomy JSON to {output_path}: {exc}"
    return None


def read_taxonomy(
    short_name: str,
    accession: str,
    metadata_dir: Path,
) -> Tuple[dict, Optional[str]]:
    """
    Read taxonomy from metadata/taxonomy_{short_name}_{accession}.json.
    If the file is absent, attempt to download it first.
    Returns (rank_dict, error_message_or_None).
    """
    tax_path = metadata_dir / f"taxonomy_{short_name}_{accession}.json"

    if not tax_path.exists():
        # Attempt to download
        tax_id = get_tax_id_from_metadata(short_name, accession, metadata_dir)
        if not tax_id:
            return {}, f"tax_id not found for {short_name}_{accession}; skipping taxonomy"
        print(
            f"  Downloading taxonomy for {short_name} (tax_id={tax_id})...",
            file=sys.stderr,
        )
        err = download_taxonomy_json(tax_id, tax_path)
        if err:
            return {}, err

    data, err = load_json_file(tax_path)
    if err or not data:
        return {}, err

    reports = data.get("reports") or []
    if not reports:
        return {}, f"No reports in taxonomy JSON: {tax_path}"
    classification = (
        (reports[0].get("taxonomy") or {}).get("classification") or {}
    )
    result = {}
    for rank in TAXONOMY_RANKS:
        entry = classification.get(rank)
        if entry and entry.get("name"):
            result[rank] = entry["name"]
    if not result:
        return {}, f"No classification data in {tax_path}"
    return result, None


# ---------------------------------------------------------------------------
# Row building
# ---------------------------------------------------------------------------

def build_org_rows(
    trio_id: str,
    manifest: dict,
    run_dir: Path,
) -> Tuple[List[OrgRow], List[str]]:
    """Build one OrgRow per organism entry in the manifest."""
    metadata_dir = run_dir / "metadata"
    noncoding = has_noncoding_data(run_dir)
    rows: List[OrgRow] = []
    warnings: List[str] = []

    for org in manifest.get("organisms", []):
        short_name = org.get("short_name", "")
        accession = org.get("accession", "")
        role = org.get("role", "")
        organism_name = org.get("raw_organism_name") or org.get("ncbi_full_name", "")

        tax, err = read_taxonomy(short_name, accession, metadata_dir)
        notes = ""
        if err:
            warnings.append(f"[{trio_id}/{short_name}] {err}")
            notes = err

        rows.append(OrgRow(
            trio_id=trio_id,
            has_noncoding=noncoding,
            organism_name=organism_name,
            accession=accession,
            role=role,
            kingdom=tax.get("kingdom", MISSING),
            phylum=tax.get("phylum", MISSING),
            klass=tax.get("class", MISSING),
            order=tax.get("order", MISSING),
            family=tax.get("family", MISSING),
            genus=tax.get("genus", MISSING),
            species=tax.get("species", MISSING),
            notes=notes,
        ))

    return rows, warnings


def collect_all_rows(input_dir: Path) -> Tuple[List[OrgRow], List[str]]:
    """Walk all trio directories and build OrgRows, accumulating warnings."""
    all_rows: List[OrgRow] = []
    all_warnings: List[str] = []

    trio_dirs = find_trio_dirs(input_dir)
    if not trio_dirs:
        all_warnings.append(f"No trio directories found under {input_dir}")
        return all_rows, all_warnings

    for trio_dir in trio_dirs:
        trio_id = trio_dir.name
        run_dir = find_latest_run_dir(trio_dir)
        if run_dir is None:
            all_warnings.append(f"[{trio_id}] No dated run directory found; skipping")
            continue

        manifest, err = load_manifest(run_dir)
        if err or manifest is None:
            all_warnings.append(f"[{trio_id}] {err}; skipping")
            continue

        rows, warnings = build_org_rows(trio_id, manifest, run_dir)
        all_rows.extend(rows)
        all_warnings.extend(warnings)

    return all_rows, all_warnings


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def _cell(text: str) -> str:
    """Escape pipe characters in a Markdown table cell."""
    return str(text).replace("|", "\\|")


def render_markdown_table(rows: List[OrgRow]) -> str:
    headers = [
        "Trio", "Role", "Organism name", "Accession",
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    ]
    sep = [":---", ":---", ":---", ":---",
           ":---", ":---", ":---", ":---", ":---", ":---", ":---"]

    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(sep) + " |",
    ]

    prev_trio_id = None
    for r in rows:
        # Show trio ID with Non-coding annotation only on the first row of each trio.
        if r.trio_id != prev_trio_id:
            ncds_label = "Non-coding: Yes" if r.has_noncoding else "Non-coding: No"
            trio_cell = f"**{_cell(r.trio_id)}** *[{ncds_label}]*"
            prev_trio_id = r.trio_id
        else:
            trio_cell = ""

        cols = [
            trio_cell,
            _cell(r.role),
            _cell(r.organism_name),
            _cell(r.accession),
            _cell(r.kingdom),
            _cell(r.phylum),
            _cell(r.klass),
            _cell(r.order),
            _cell(r.family),
            _cell(r.genus),
            _cell(r.species),
        ]
        lines.append("| " + " | ".join(cols) + " |")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize all trios under a lineage directory into a Markdown table."
    )
    parser.add_argument(
        "input_dir",
        help="Lineage results directory (e.g., results/fungi/).",
    )
    parser.add_argument(
        "--output",
        help=(
            "Output Markdown file path. "
            "Defaults to <input_dir>/trio_summary.md. "
            "Use '-' for stdout."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Determine output path
    if args.output and args.output != "-":
        output_path: Optional[Path] = Path(args.output)
    elif args.output == "-":
        output_path = None  # stdout
    else:
        output_path = input_dir / "trio_summary.md"

    rows, warnings = collect_all_rows(input_dir)

    for w in warnings:
        print(f"Warning: {w}", file=sys.stderr)

    if not rows:
        print("No rows to output.", file=sys.stderr)
        sys.exit(1)

    table = render_markdown_table(rows)

    if output_path is None:
        print(table)
    else:
        try:
            output_path.write_text(table, encoding="utf-8")
            print(f"Written: {output_path}", file=sys.stderr)
        except OSError as exc:
            print(f"Error writing {output_path}: {exc}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
