#!/usr/bin/env python3

"""
Fetch NCBI assembly metadata for a list of genome accessions, one row per accession.

trio_selection.R uses this to prune a phylogenetic tree down to one assembly per
species before any genome is downloaded, so the expensive stages only ever see
one representative per species.

    python3 src/select/fetch_assembly_metadata.py --accessions accs.txt --out meta.tsv

Accessions may come from a file (one per line) or from --accession, repeated.
The output TSV carries the NCBI status, taxonomy, assembly structure, ANI,
CheckM, and annotation BUSCO fields needed to audit and rank Stage-0 candidates.
Accessions that NCBI no longer returns are reported on stderr and omitted from
the table -- callers must treat a missing row as "drop this leaf", not as an
error.
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
    "current_accession",
    "organism_name",
    "organism_tax_id",
    "species",
    "genus",
    "source_database",
    "refseq_category",
    "assembly_level",
    "assembly_type",
    "release_date",
    "sequencing_tech",
    "assembly_method",
    "is_atypical",
    "atypical_warnings",
    "ani_check_status",
    "checkm_completeness",
    "checkm_contamination",
    "busco_lineage",
    "busco_version",
    "busco_complete",
    "busco_duplicated",
    "busco_fragmented",
    "busco_missing",
    "total_sequence_length",
    "total_ungapped_length",
    "number_of_contigs",
    "contig_n50",
    "number_of_scaffolds",
    "scaffold_n50",
    "has_annotation",
    "annotation_provider",
    "annotation_release_date",
    "assembly_status",
    "paired_accession",
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


def build_taxon_payload(taxons: List[str], page_token: str = "") -> dict:
    """Build the reproducible NCBI request used for Stage-0 candidate discovery."""
    payload = {
        "taxons": taxons,
        "filters": {"assembly_version": "current", "mag": "exclude"},
        "page_size": PAGE_SIZE,
        "tax_exact_match": True,
    }
    if page_token:
        payload["page_token"] = page_token
    return payload


def fetch_reports_by_taxon(taxons: List[str]) -> List[dict]:
    """All current non-MAG assemblies for each taxon (many reports per taxon)."""
    reports: List[dict] = []
    page_token = ""
    while True:
        payload = build_taxon_payload(taxons, page_token)
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
    species = " ".join(tokens)
    return {"genus": genus, "species": species}


def stringify(value) -> str:
    """Convert a JSON scalar to TSV text without turning numeric zero into missing."""
    return "" if value is None else str(value)


def report_to_row(report: dict) -> Dict[str, str]:
    organism = report.get("organism") or {}
    assembly_info = report.get("assembly_info") or {}
    assembly_stats = report.get("assembly_stats") or {}
    atypical = assembly_info.get("atypical") or {}
    ani = report.get("average_nucleotide_identity") or {}
    checkm = report.get("checkm_info") or {}
    annotation = report.get("annotation_info") or {}
    busco = annotation.get("busco") or {}

    organism_name = organism.get("organism_name") or ""
    names = split_organism_name(organism_name)

    return {
        "accession": report.get("accession") or "",
        "current_accession": report.get("current_accession") or "",
        "organism_name": organism_name,
        "organism_tax_id": stringify(organism.get("tax_id")),
        "species": names["species"],
        "genus": names["genus"],
        "source_database": report.get("source_database") or "",
        "refseq_category": assembly_info.get("refseq_category") or "",
        "assembly_level": assembly_info.get("assembly_level") or "",
        "assembly_type": assembly_info.get("assembly_type") or "",
        "release_date": assembly_info.get("release_date") or "",
        "sequencing_tech": assembly_info.get("sequencing_tech") or "",
        "assembly_method": assembly_info.get("assembly_method") or "",
        "is_atypical": "true" if atypical.get("is_atypical") is True else "false",
        "atypical_warnings": "; ".join(atypical.get("warnings") or []),
        "ani_check_status": ani.get("taxonomy_check_status") or "",
        "checkm_completeness": stringify(checkm.get("completeness")),
        "checkm_contamination": stringify(checkm.get("contamination")),
        "busco_lineage": busco.get("busco_lineage") or "",
        "busco_version": busco.get("busco_ver") or "",
        "busco_complete": stringify(busco.get("complete")),
        "busco_duplicated": stringify(busco.get("duplicated")),
        "busco_fragmented": stringify(busco.get("fragmented")),
        "busco_missing": stringify(busco.get("missing")),
        "total_sequence_length": stringify(assembly_stats.get("total_sequence_length")),
        "total_ungapped_length": stringify(assembly_stats.get("total_ungapped_length")),
        "number_of_contigs": stringify(assembly_stats.get("number_of_contigs")),
        "contig_n50": stringify(assembly_stats.get("contig_n50")),
        "number_of_scaffolds": stringify(assembly_stats.get("number_of_scaffolds")),
        "scaffold_n50": stringify(assembly_stats.get("scaffold_n50")),
        "has_annotation": "true" if annotation else "false",
        "annotation_provider": annotation.get("provider") or "",
        "annotation_release_date": annotation.get("release_date") or "",
        "assembly_status": assembly_info.get("assembly_status") or "",
        "paired_accession": report.get("paired_accession") or "",
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
        help="File with one species name or NCBI Taxonomy ID per line; fetch all current "
        "non-MAG assemblies for each taxon (instead of specific accessions).",
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
