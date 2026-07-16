#!/usr/bin/env python3

"""
Fetch NCBI assembly metadata for a list of genome accessions, one row per accession.

trio_selection.R uses this to prune a phylogenetic tree down to one assembly per
species before any genome is downloaded, so the expensive stages only ever see
one representative per species.

    python3 src/select/fetch_assembly_metadata.py --accessions accs.txt --out meta.tsv

Accessions may come from a file (one per line) or from --accession, repeated.
The output TSV carries accession, organism_name, species, genus, refseq_category,
assembly_level, contig_n50 and assembly_status.  Accessions that NCBI no longer
returns are reported on stderr and omitted from the table -- callers must treat a
missing row as "drop this leaf", not as an error.
"""

import argparse
import csv
import json
import sys
import urllib.error
import urllib.request
from typing import Dict, Iterable, List

DATASET_REPORT_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"

# The v2 endpoint defaults to page_size=20 and truncates *silently*, so the page
# size must be set explicitly and the response paged until it is exhausted.
PAGE_SIZE = 1000

COLUMNS = (
    "accession",
    "organism_name",
    "species",
    "genus",
    "refseq_category",
    "assembly_level",
    "contig_n50",
    "total_ungapped_length",
    "has_annotation",
    "assembly_status",
)


def post_json(url: str, payload: dict) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_reports(accessions: List[str]) -> List[dict]:
    reports: List[dict] = []
    page_token = ""
    while True:
        payload = {"accessions": accessions, "page_size": PAGE_SIZE}
        if page_token:
            payload["page_token"] = page_token
        response = post_json(DATASET_REPORT_URL, payload)
        reports.extend(response.get("reports", []))
        page_token = response.get("next_page_token", "")
        if not page_token:
            break
    return reports


def fetch_reports_by_taxon(taxons: List[str]) -> List[dict]:
    """All current, non-atypical assemblies for each taxon (many reports per taxon)."""
    reports: List[dict] = []
    page_token = ""
    while True:
        payload = {
            "taxons": taxons,
            "filters": {"assembly_version": "current", "exclude_atypical": True},
            "page_size": PAGE_SIZE,
        }
        if page_token:
            payload["page_token"] = page_token
        response = post_json(DATASET_REPORT_URL, payload)
        reports.extend(response.get("reports", []))
        page_token = response.get("next_page_token", "")
        if not page_token:
            break
    return reports


def read_lines(path: str) -> List[str]:
    """Read non-blank, non-comment lines from a file, de-duplicated, order-preserving."""
    with open(path, encoding="utf-8") as handle:
        items = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    return list(dict.fromkeys(items))


def split_organism_name(organism_name: str) -> Dict[str, str]:
    tokens = organism_name.split()
    genus = tokens[0] if tokens else ""
    # Subspecies and strain suffixes are dropped: two assemblies of the same
    # species must collapse to the same key even when one is named to subspecies.
    species = " ".join(tokens[:2]) if len(tokens) >= 2 else genus
    return {"genus": genus, "species": species}


def report_to_row(report: dict) -> Dict[str, str]:
    organism = report.get("organism") or {}
    assembly_info = report.get("assembly_info") or {}
    assembly_stats = report.get("assembly_stats") or {}

    organism_name = organism.get("organism_name") or ""
    names = split_organism_name(organism_name)

    return {
        "accession": report.get("accession") or "",
        "organism_name": organism_name,
        "species": names["species"],
        "genus": names["genus"],
        "refseq_category": assembly_info.get("refseq_category") or "",
        "assembly_level": assembly_info.get("assembly_level") or "",
        "contig_n50": str(assembly_stats.get("contig_n50") or ""),
        "total_ungapped_length": str(assembly_stats.get("total_ungapped_length") or ""),
        "has_annotation": "true" if report.get("annotation_info") else "false",
        "assembly_status": assembly_info.get("assembly_status") or "",
    }


def read_accessions(args: argparse.Namespace) -> List[str]:
    accessions: List[str] = list(args.accession or [])
    if args.accessions:
        with open(args.accessions, encoding="utf-8") as handle:
            accessions.extend(
                line.strip() for line in handle if line.strip() and not line.startswith("#")
            )
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(accessions))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch NCBI assembly metadata for genome accessions as a TSV."
    )
    parser.add_argument("--accessions", help="File with one accession per line.")
    parser.add_argument(
        "--accession", action="append", help="Accession to fetch; may be repeated."
    )
    parser.add_argument(
        "--taxons",
        help="File with one species name per line; fetch ALL current, non-atypical "
        "assemblies for each species (instead of specific accessions).",
    )
    parser.add_argument("--out", help="Output TSV path. Defaults to stdout.")
    return parser.parse_args()


def write_rows(rows: Iterable[Dict[str, str]], out_path: str) -> None:
    handle = open(out_path, "w", newline="", encoding="utf-8") if out_path else sys.stdout
    try:
        writer = csv.DictWriter(
            handle, fieldnames=list(COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    finally:
        if out_path:
            handle.close()


def main() -> int:
    args = parse_args()

    if args.taxons:
        taxons = read_lines(args.taxons)
        if not taxons:
            print("error: no taxons given in --taxons file.", file=sys.stderr)
            return 1
        try:
            reports = fetch_reports_by_taxon(taxons)
        except urllib.error.URLError as exc:
            print(f"error: NCBI request failed: {exc}", file=sys.stderr)
            return 1
        rows = [report_to_row(report) for report in reports]
        present = {row["species"] for row in rows if row["assembly_status"] == "current"}
        for taxon in taxons:
            if taxon not in present:
                print(
                    f"warning: NCBI returned no current assembly for taxon '{taxon}'; "
                    "species dropped.",
                    file=sys.stderr,
                )
        write_rows(rows, args.out)
        print(
            f"{len(rows)} assemblies for {len(present)}/{len(taxons)} taxa.",
            file=sys.stderr,
        )
        return 0

    accessions = read_accessions(args)
    if not accessions:
        print("error: no accessions given (use --accessions or --accession).", file=sys.stderr)
        return 1

    try:
        reports = fetch_reports(accessions)
    except urllib.error.URLError as exc:
        print(f"error: NCBI request failed: {exc}", file=sys.stderr)
        return 1

    rows = [report_to_row(report) for report in reports]
    found = {row["accession"] for row in rows}

    missing = [accession for accession in accessions if accession not in found]
    for accession in missing:
        print(f"warning: NCBI returned no report for {accession}; dropping it.", file=sys.stderr)

    # Emit in the order requested so the caller can join positionally if it wants.
    by_accession = {row["accession"]: row for row in rows}
    ordered = [by_accession[a] for a in accessions if a in by_accession]

    write_rows(ordered, args.out)
    print(
        f"{len(ordered)}/{len(accessions)} accessions resolved.", file=sys.stderr
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
