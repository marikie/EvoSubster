import contextlib
import io
import os
import sys
import time
import unittest
import urllib.error
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "select"))

import fetch_assembly_metadata as f


LEGACY_COLUMNS = (
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
)


class ReportToRowTest(unittest.TestCase):
    def test_normalizes_complete_organism_name_whitespace(self):
        self.assertEqual(
            f.normalize_organism_name("  Genus\talpha  strain-X "),
            "Genus alpha strain-X",
        )

    def test_extracts_quality_fields_for_stage0_ranking(self):
        report = {
            "accession": "GCF_000000001.1",
            "current_accession": "GCF_000000001.1",
            "source_database": "SOURCE_DATABASE_REFSEQ",
            "organism": {"tax_id": 1234, "organism_name": "Genus alpha strain-X"},
            "assembly_info": {
                "refseq_category": "reference genome",
                "assembly_level": "Chromosome",
                "assembly_status": "current",
                "assembly_type": "haploid",
                "release_date": "2025-01-02",
                "sequencing_tech": "PacBio; Illumina",
                "assembly_method": "Assembler v1",
                "atypical": {"is_atypical": False, "warnings": []},
            },
            "assembly_stats": {
                "contig_n50": 924431,
                "number_of_contigs": 18,
                "number_of_scaffolds": 12,
                "scaffold_n50": 2000000,
                "total_sequence_length": "12500000",
                "total_ungapped_length": "12071326",
            },
            "average_nucleotide_identity": {"taxonomy_check_status": "OK"},
            "checkm_info": {"completeness": 99.4, "contamination": 0.2},
            "annotation_info": {
                "name": "GCF annotation",
                "provider": "NCBI RefSeq",
                "release_date": "2025-03-04",
                "busco": {
                    "busco_lineage": "example_odb10",
                    "busco_ver": "5.7.1",
                    "complete": 0.98,
                    "duplicated": 0.01,
                    "fragmented": 0.005,
                    "missing": 0.015,
                },
            },
        }
        row = f.report_to_row(report)
        self.assertEqual(row["species"], "Genus alpha strain-X")
        self.assertEqual(row["organism_tax_id"], "1234")
        self.assertEqual(row["source_database"], "SOURCE_DATABASE_REFSEQ")
        self.assertEqual(row["contig_n50"], "924431")
        self.assertEqual(row["total_ungapped_length"], "12071326")
        self.assertEqual(row["ani_check_status"], "OK")
        self.assertEqual(row["checkm_completeness"], "99.4")
        self.assertEqual(row["busco_lineage"], "example_odb10")
        self.assertEqual(row["busco_complete"], "0.98")
        self.assertEqual(row["is_atypical"], "false")
        self.assertEqual(row["has_annotation"], "true")
        self.assertEqual(row["assembly_status"], "current")

    def test_missing_annotation_and_stats(self):
        report = {
            "accession": "GCA_000000002.1",
            "organism": {"organism_name": "Genus beta"},
            "assembly_info": {"assembly_level": "Contig", "assembly_status": "current"},
            "assembly_stats": {"contig_n50": 5000},
        }
        row = f.report_to_row(report)
        self.assertEqual(row["has_annotation"], "false")
        self.assertEqual(row["total_ungapped_length"], "")
        self.assertEqual(row["paired_accession"], "")

    def test_records_paired_accession_from_dataset_report(self):
        report = {
            "accession": "GCF_901000725.3",
            "paired_accession": "GCA_901000725.3",
            "organism": {"organism_name": "Genus paired"},
        }
        row = f.report_to_row(report)
        self.assertEqual(row["paired_accession"], "GCA_901000725.3")
        self.assertEqual(tuple(row.keys()), f.COLUMNS)

    def test_appends_paired_accession_without_moving_public_columns(self):
        self.assertEqual(f.COLUMNS[:-1], LEGACY_COLUMNS)
        self.assertEqual(f.COLUMNS[-1], "paired_accession")

    def test_preserves_atypical_warning_for_local_audit(self):
        report = {
            "accession": "GCA_000000003.1",
            "organism": {"tax_id": 4321, "organism_name": "Genus delta"},
            "assembly_info": {
                "assembly_status": "current",
                "atypical": {
                    "is_atypical": True,
                    "warnings": ["contaminated", "genome length too large"],
                },
            },
        }
        row = f.report_to_row(report)
        self.assertEqual(row["is_atypical"], "true")
        self.assertEqual(row["atypical_warnings"], "contaminated; genome length too large")

    def test_columns_match_dictwriter_fieldnames(self):
        # Every key report_to_row emits must be a declared column, and vice versa.
        report = {"accession": "GCA_1", "organism": {"organism_name": "Genus gamma"},
                  "assembly_info": {}, "assembly_stats": {}}
        self.assertEqual(tuple(f.report_to_row(report).keys()), f.COLUMNS)


class TaxonRequestTest(unittest.TestCase):
    def test_retries_a_transient_bad_gateway_response(self):
        transient_error = urllib.error.HTTPError(
            f.DATASET_REPORT_URL,
            502,
            "Bad Gateway",
            {},
            io.BytesIO(b"temporary upstream failure"),
        )
        success_response = io.BytesIO(b'{"reports": []}')

        with mock.patch.object(
            f.urllib.request,
            "urlopen",
            side_effect=[transient_error, success_response],
        ) as urlopen, mock.patch.object(time, "sleep"):
            with contextlib.redirect_stderr(io.StringIO()):
                try:
                    response = f.post_json(f.DATASET_REPORT_URL, {"taxons": ["9606"]})
                except urllib.error.HTTPError as exc:
                    self.fail(
                        f"post_json did not retry the transient HTTP {exc.code} response"
                    )

        self.assertEqual(response, {"reports": []})
        self.assertEqual(urlopen.call_count, 2)

    def test_requests_current_non_mag_candidates_for_local_quality_audit(self):
        payload = f.build_taxon_payload(["9606", "562"])
        self.assertEqual(
            payload,
            {
                "taxons": ["9606", "562"],
                "filters": {"assembly_version": "current", "mag": "exclude"},
                "page_size": f.PAGE_SIZE,
                "tax_exact_match": True,
            },
        )


if __name__ == "__main__":
    unittest.main()
