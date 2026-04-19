#!/usr/bin/env python3

"""
Collect genome sizes and pairwise identity values for each benchmark trio.

Reads per-lineage summary JSONs to obtain pairwise identity (idt_12/13/23) and
species accessions, then looks up total assembly sizes from NCBI metadata JSON
files. Outputs a TSV keyed by trio_label compatible with bench_summary.tsv.

Usage:
    python collect_bench_metadata.py <results_dir> <genomes_dir> [--output PATH]

    results_dir  - directory containing <lineage>/<lineage>_summary.json files
    genomes_dir  - root of NCBI genome JSON files (searched recursively for
                   {accession}.json)
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Optional


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", help="Directory containing lineage summary JSONs")
    parser.add_argument("genomes_dir", help="Root directory of NCBI genome JSON files")
    parser.add_argument(
        "--output",
        default=None,
        help="Output TSV path (default: <results_dir>/benchmark/bench_metadata.tsv)",
    )
    return parser.parse_args()


def build_accession_index(genomes_dir: Path) -> Dict[str, Path]:
    """Return a mapping from accession ID to its JSON file path."""
    index: Dict[str, Path] = {}
    for p in genomes_dir.rglob("*.json"):
        # Match files named exactly like GCA_XXXXXXX.N.json or GCF_XXXXXXX.N.json
        if p.stem.startswith("GCA_") or p.stem.startswith("GCF_"):
            index[p.stem] = p
    return index


def get_genome_size(accession: str, acc_index: Dict[str, Path]) -> Optional[int]:
    """Return total_sequence_length (bp) from the NCBI metadata JSON, or None."""
    p = acc_index.get(accession)
    if p is None:
        return None
    try:
        data = json.loads(p.read_text())
        reports = data.get("reports", [])
        if not reports:
            return None
        return reports[0].get("assembly_stats", {}).get("total_sequence_length")
    except Exception:
        return None


def short_name_to_label(short_name: str) -> str:
    """Convert a summary-JSON short name (e.g. 'Agabis1') to a bench label (e.g. 'agaBis')."""
    return short_name[:3].lower() + short_name[3:6].capitalize()


def parse_identity(value: str) -> float:
    """Parse an identity string like '88.7999 %' into a float."""
    return float(str(value).replace("%", "").strip())


HEADER = [
    "trio_label",
    "lineage",
    "size_org1_bp",
    "size_org2_bp",
    "size_org3_bp",
    "idt_12",
    "idt_13",
    "idt_23",
]


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    genomes_dir = Path(args.genomes_dir)

    if args.output:
        out_path = Path(args.output)
    else:
        out_path = results_dir / "benchmark" / "bench_metadata.tsv"

    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Building accession index from {genomes_dir} ...", file=sys.stderr)
    acc_index = build_accession_index(genomes_dir)
    print(f"  Found {len(acc_index)} genome JSON files", file=sys.stderr)

    # Collect lineage summary JSONs only (parent dir is all-lowercase, e.g. fungi, cnidaria).
    # Exclude dataset-level summary JSONs whose parent dir is a dataset name (contains capitals).
    summary_files = sorted(
        p for p in results_dir.glob("*/*_summary.json") if p.parent.name.islower()
    )
    if not summary_files:
        print(f"ERROR: no *_summary.json files found under {results_dir}", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing_size = []

    for sf in summary_files:
        lineage = sf.parent.name
        try:
            data = json.loads(sf.read_text())
        except Exception as e:
            print(f"WARNING: could not read {sf}: {e}", file=sys.stderr)
            continue

        for entry in data.get("datasets", []):
            if entry.get("filter_status") == "filtered":
                continue

            sn1 = entry["species1"]["metadata"]["short_name"]
            sn2 = entry["species2"]["metadata"]["short_name"]
            sn3 = entry["species3"]["metadata"]["short_name"]
            trio_label = "_".join(short_name_to_label(s) for s in [sn1, sn2, sn3])

            acc1 = entry["species1"].get("accession", "")
            acc2 = entry["species2"].get("accession", "")
            acc3 = entry["species3"].get("accession", "")

            size1 = get_genome_size(acc1, acc_index)
            size2 = get_genome_size(acc2, acc_index)
            size3 = get_genome_size(acc3, acc_index)

            for label, acc, sz in [(sn1, acc1, size1), (sn2, acc2, size2), (sn3, acc3, size3)]:
                if sz is None:
                    missing_size.append(f"{trio_label}/{label} ({acc})")

            idt_12 = parse_identity(entry.get("idt_12", "NA"))
            idt_13 = parse_identity(entry.get("idt_13", "NA"))
            idt_23 = parse_identity(entry.get("idt_23", "NA"))

            rows.append(
                {
                    "trio_label": trio_label,
                    "lineage": lineage,
                    "size_org1_bp": size1 if size1 is not None else "NA",
                    "size_org2_bp": size2 if size2 is not None else "NA",
                    "size_org3_bp": size3 if size3 is not None else "NA",
                    "idt_12": f"{idt_12:.4f}",
                    "idt_13": f"{idt_13:.4f}",
                    "idt_23": f"{idt_23:.4f}",
                }
            )

    if missing_size:
        print("WARNING: genome size not found for:", file=sys.stderr)
        for m in missing_size:
            print(f"  {m}", file=sys.stderr)

    rows.sort(key=lambda r: (r["lineage"], r["trio_label"]))

    with out_path.open("w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[col]) for col in HEADER) + "\n")

    print(f"Wrote {len(rows)} rows to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
