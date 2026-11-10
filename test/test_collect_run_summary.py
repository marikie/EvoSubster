import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "report"))

import collect_run_summary as script
import report_contract


SLOTS = (
    ("org1", "Out1", "GCF_1", "Pan alpha"),
    ("org2", "In2", "GCA_2", "Homo beta"),
    ("org3", "In3", "GCA_3", "Homo gamma"),
)


def write_file(path, content=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def create_run(dataset_dir, layout="current", date="20260101", stale_metadata=False):
    run = dataset_dir / date
    organisms = []
    artifacts = []
    def artifact(kind, slot, relative, content):
        path = run / relative
        write_file(path, content)
        entry = dict(kind=kind, slot=slot, path=relative, sha256=script.runs.sha256(path))
        artifacts.append(entry)
        return entry
    for slot, name, accession, organism in SLOTS:
        meta = f"metadata/{slot}.json"
        tax = f"metadata/{slot}-taxonomy.json"
        write_file(run / meta, "{}")
        write_file(run / tax, json.dumps({'reports':[{'taxonomy':{'classification':{
            'genus':{'name': organism.split()[0]}, 'species':{'name': organism}}}}]}))
        organisms.append(dict(slot=slot, role='outgroup' if slot == 'org1' else 'ingroup',
            short_name=name, accession=accession, raw_organism_name=organism,
            fasta_sha256='fixture', gff_sha256=None, metadata_json=meta, taxonomy_json=tax,
            metadata_json_sha256=script.runs.sha256(run / meta),
            taxonomy_json_sha256=script.runs.sha256(run / tax)))
        artifact('gc', slot, f"arbitrary/{slot}-gc.txt", "Total GC content: 40%\n")
        artifact('coding', slot, f"arbitrary/{slot}-coding.txt", "Coding region percentage: 25.00%\n")
        if slot != 'org1':
            artifact('single', slot, f"tables/{slot}.tsv", "context\tcount\n")
            artifact('dinuc', slot, f"tables/{slot}-dinuc.tsv", "mutType\tmutNum\ttotalRootNum\n")
            for kind in ('norm','sbst','ori','logratio','dinuc_tsv'):
                artifact('pdf_'+kind, slot, f"plots/{slot}-{kind}.pdf", "%PDF-fixture\n%%EOF")
    artifact('sbst_ratio', 'all', 'ratios.txt', 'org2: 1.2\norg3: 1.3\n')
    for pair, value in (('org1-org2',90),('org1-org3',91),('org2-org3',97)):
        content = '# substitution percent identity: '+str(value)+'\n# substitution percent identity: 99\n'
        content += ''.join('#last -'+option+' 1\n' for option in ('t','a','A','b','B','S'))
        content += '# score matrix (query letters = columns, reference letters = rows):\n'
        content += ''.join(letter+' 1 0 0 0\n' for letter in 'ACGT')
        entry = artifact('train', pair, f"training/{pair}.txt", content)
        entry.update(identity_percent=value, identity_source="last-train-penultimate", reference_slot=pair.split('-')[0], query_slot=pair.split('-')[1])
    payload = dict(schema_version=2, trio_key='fixture-trio', run_key='fixture-run',
                   created_at=date, date=date, options={}, status='complete',
                   organisms=organisms, artifacts=artifacts)
    if layout == 'legacy': payload.pop('schema_version')
    write_file(run / script.runs.MANIFEST, json.dumps(payload))
    return run


class BuildSummaryManifestTest(unittest.TestCase):
    def test_reads_explicit_artifacts_and_final_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = create_run(Path(tmp)/'arbitrary')
            summary, _ = script.build_summary(Path(tmp),80)
            entry = summary['datasets'][0]
            self.assertEqual(entry['dataset'],'fixture-trio')
            self.assertEqual(entry['identity_source'],'last-train-penultimate')
            self.assertEqual(entry['idt_12'],'90 %')
            self.assertEqual(entry['species2']['tsv_file'],str(run/'tables/org2.tsv'))
            self.assertEqual(entry['species1']['coding_region_percentage'],['Coding region percentage: 25.00%'])
            self.assertEqual(report_contract.validate_summary(summary)[1],[])

    def test_trio_table_uses_recorded_identity(self):
        import summarize_trios
        with tempfile.TemporaryDirectory() as tmp:
            create_run(Path(tmp)/'arbitrary-folder')
            rows,warnings=summarize_trios.collect_all_rows(Path(tmp))
            self.assertEqual(len(rows),3)
            self.assertEqual({row.trio_id for row in rows},{'fixture-trio'})

    def test_directory_rename_preserves_logical_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = create_run(root/'old')
            first, _ = script.build_summary(root,80)
            run.rename(root/'renamed')
            second, _ = script.build_summary(root,80)
            for field in ('dataset','latest_run','idt_12','idt_13','idt_23'):
                self.assertEqual(first['datasets'][0][field],second['datasets'][0][field])

    def test_legacy_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            create_run(Path(tmp)/'old', 'legacy')
            with self.assertRaises(ValueError): script.build_summary(Path(tmp),80)

    def test_modified_artifact_is_rejected_without_filename_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = create_run(Path(tmp)/'new')
            (run/'tables/org2.tsv').write_text('modified')
            with self.assertRaises(ValueError): script.build_summary(Path(tmp),80)

    def test_incomplete_run_is_not_a_full_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = create_run(Path(tmp)/'new')
            payload = script.runs.load(run)
            payload['status']='identity_complete'
            script.runs.atomic_json(run/script.runs.MANIFEST,payload)
            summary,_=script.build_summary(Path(tmp),80)
            self.assertEqual(summary['datasets'],[])



class ReportContractTest(unittest.TestCase):
    def test_derives_default_and_filtered_summary_paths(self):
        root = Path("/tmp/results/fungi")
        summary, filtered = report_contract.summary_paths(root)
        self.assertEqual(summary, root / "fungi_summary.json")
        self.assertEqual(filtered, root / "fungi_summary_filtered.json")

    def test_report_extension_matches_output_format(self):
        summary = Path("/tmp/fungi_summary.json")
        self.assertEqual(
            report_contract.report_output_path(summary, "word_document"),
            Path("/tmp/fungi_summary.docx"),
        )
        self.assertEqual(
            report_contract.report_output_path(summary, "html_document"),
            Path("/tmp/fungi_summary.html"),
        )
        self.assertEqual(
            report_contract.report_output_path(summary, "pdf_document"),
            Path("/tmp/fungi_summary.pdf"),
        )

    def test_rejects_unsupported_output_format(self):
        with self.assertRaisesRegex(ValueError, "unsupported output format"):
            report_contract.report_output_path(
                Path("/tmp/fungi_summary.json"), "github_document"
            )

    def test_validates_required_dataset_fields(self):
        dataset = {"dataset": "TrioA"}
        names, errors = report_contract.validate_summary({"datasets": [dataset]})
        self.assertEqual(names, ["TrioA"])
        self.assertEqual(len(errors), 1)
        self.assertIn("species1", errors[0])
        self.assertIn("idt_23", errors[0])

    def test_rejects_invalid_species_and_pdf_types(self):
        dataset = {
            "dataset": "TrioA",
            "species1": {"metadata": "not-an-object"},
            "species2": "not-an-object",
            "species3": {
                "taxonomy": "not-an-object",
                "gc_content": "not-an-array",
                "pdfs": {"norm": "not-an-array"},
            },
            "idt_12": "90 %",
            "idt_13": 91,
            "idt_23": "97 %",
        }

        _, errors = report_contract.validate_summary({"datasets": [dataset]})

        self.assertTrue(any("species2 must be an object" in error for error in errors))
        self.assertTrue(any("metadata must be an object" in error for error in errors))
        self.assertTrue(any("taxonomy must be an object" in error for error in errors))
        self.assertTrue(any("gc_content must be an array" in error for error in errors))
        self.assertTrue(any("pdfs.norm must be an array" in error for error in errors))
        self.assertTrue(any("idt_13 must be a string" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
