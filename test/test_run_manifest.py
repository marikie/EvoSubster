import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
import run_manifest as manifest


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / 'input.fa'
        self.fasta.write_text('>chr1\nACGT\n')

    def test_bundle_detects_changed_fasta_and_gff(self):
        gff = self.root / 'genomic.gff'
        gff.write_text('chr1\tx\tCDS\t1\t3\t.\t+\t0\tID=x\n')
        manifest.bundle_save(self.root, 'GCA_1.1', self.fasta)
        self.assertEqual(manifest.bundle_gff(self.fasta), str(gff))
        gff.write_text('altered')
        with self.assertRaises(ValueError):
            manifest.bundle_gff(self.fasta)

    def test_bundle_missing_annotation_is_explicit(self):
        manifest.bundle_save(self.root, 'GCA_1.1', self.fasta)
        (self.root / 'genomic.gff').write_text('stale')
        self.assertEqual(manifest.bundle_gff(self.fasta), 'NO_GFF_FILE')

    def test_bundle_rejects_wrong_sequence_ids(self):
        (self.root / 'genomic.gff').write_text('other\tx\tCDS\t1\t3\t.\t+\t0\tID=x\n')
        with self.assertRaises(ValueError):
            manifest.bundle_save(self.root, 'GCA_1.1', self.fasta)

    def test_training_key_tracks_content_direction_options(self):
        query = self.root / 'query.fa'
        query.write_text('>chr2\nAGGT\n')
        with patch.object(manifest.subprocess, 'check_output', return_value='lastal 1648'):
            first = manifest.training_directory(self.root, self.fasta, query, 'a', 'b')
            reverse = manifest.training_directory(self.root, query, self.fasta, 'b', 'a')
            self.assertNotEqual(first, reverse)
            self.fasta.write_text('>chr1\nAAAA\n')
            self.assertNotEqual(first, manifest.training_directory(self.root, self.fasta, query, 'a', 'b'))

    def test_checkpoint_recovery_and_relocation(self):
        run = self.root / 'run'
        value = dict(schema_version=2, run_key='x', trio_key='t', options={}, created_at='now',
                     organisms=[dict(slot='org'+str(i), fasta_sha256='x', gff_sha256=None) for i in (1,2,3)],
                     artifacts=[], status='pending')
        manifest.prepare(run, value)
        counts = run / 'statistics' / 'counts.tsv'
        counts.parent.mkdir()
        counts.write_text('validated')
        manifest.checkpoint(run, 'counts_complete', ['single:org2:statistics/counts.tsv'])
        partial = run / 'intermediateFiles' / 'db' / 'partial'
        partial.parent.mkdir(parents=True)
        partial.write_text('incomplete')
        manifest.prepare(run, value)
        self.assertTrue(counts.exists())
        self.assertFalse(partial.parent.exists())
        relocated = self.root / 'arbitrary'
        run.rename(relocated)
        self.assertEqual(manifest.load(relocated)['artifacts'][0]['path'], 'statistics/counts.tsv')

    def test_modified_metadata_cannot_be_blessed_on_resume(self):
        value = dict(schema_version=2, run_key='x', trio_key='t', options={}, created_at='now',
                     organisms=[dict(slot='org'+str(i), fasta_sha256='x', gff_sha256=None) for i in (1,2,3)],
                     artifacts=[], status='pending')
        run=self.root/'run'
        manifest.prepare(run,value)
        (run/'metadata/org1_assembly.json').write_text('tampered')
        with self.assertRaisesRegex(ValueError,'Metadata checksum mismatch'):
            manifest.prepare(run,value)
        self.assertEqual((run/'metadata/org1_assembly.json').read_text(),'tampered')

    def test_complete_manifest_requires_count_channels(self):
        from test_collect_run_summary import create_run
        run=create_run(self.root/'counts')
        value=manifest.load(run)
        value['artifacts']=[a for a in value['artifacts'] if not (a['kind']=='dinuc' and a['slot']=='org2')]
        manifest.atomic_json(run/manifest.MANIFEST,value)
        with self.assertRaisesRegex(ValueError,'Missing completed artifacts'):
            manifest.validated(run,complete=True)

    def test_duplicate_logical_artifact_is_rejected(self):
        from test_collect_run_summary import create_run
        run=create_run(self.root/'duplicates')
        value=manifest.load(run)
        duplicate=dict(next(a for a in value['artifacts'] if a['kind']=='single'))
        duplicate['path']='other.tsv'
        value['artifacts'].append(duplicate)
        manifest.atomic_json(run/manifest.MANIFEST,value)
        with self.assertRaisesRegex(ValueError,'Duplicate logical artifact'):
            manifest.load(run)

    def test_database_inventory_detects_partial_damage(self):
        db=self.root/'db'
        db.mkdir()
        (db/'one').write_text('one'); (db/'two').write_text('two')
        manifest.seal_directory(db)
        manifest.seal_directory(db,verify=True)
        (db/'one').unlink()
        with self.assertRaises(ValueError): manifest.seal_directory(db,verify=True)

    def test_unmanifested_output_rejected(self):
        run = self.root / 'run'
        run.mkdir()
        (run / 'old-data').write_text('old')
        with self.assertRaises(ValueError):
            manifest.prepare(run, dict(run_key='x'))
        self.assertEqual((run / 'old-data').read_text(), 'old')
