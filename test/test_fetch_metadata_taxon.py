import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "select"))

import fetch_assembly_metadata as f


class ReportToRowTest(unittest.TestCase):
    def test_extracts_new_columns_and_annotation(self):
        report = {
            "accession": "GCF_000000001.1",
            "organism": {"organism_name": "Genus alpha strain-X"},
            "assembly_info": {
                "refseq_category": "reference genome",
                "assembly_level": "Chromosome",
                "assembly_status": "current",
            },
            "assembly_stats": {
                "contig_n50": 924431,
                "total_ungapped_length": "12071326",
                "total_number_of_chromosomes": 16,
            },
            "annotation_info": {"name": "GCF annotation"},
        }
        row = f.report_to_row(report)
        self.assertEqual(row["species"], "Genus alpha")
        self.assertEqual(row["contig_n50"], "924431")
        self.assertEqual(row["total_ungapped_length"], "12071326")
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

    def test_columns_match_dictwriter_fieldnames(self):
        # Every key report_to_row emits must be a declared column, and vice versa.
        report = {"accession": "GCA_1", "organism": {"organism_name": "Genus gamma"},
                  "assembly_info": {}, "assembly_stats": {}}
        self.assertEqual(set(f.report_to_row(report).keys()), set(f.COLUMNS))


if __name__ == "__main__":
    unittest.main()
