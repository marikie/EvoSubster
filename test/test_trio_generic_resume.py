"""Local --trios acceptance coverage for generic stage reuse and read-only previews."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import unittest

import test_coding_region_pipeline as coding_fixture


@unittest.skipUnless(shutil.which('yq') and shutil.which('jq'), 'pipeline requires yq and jq')
class TrioGenericResumeTests(unittest.TestCase):
    def setUp(self):
        coding_fixture.CodingRegionPipelineTest.setUp(self)
        # Keep real entrypoints, helpers, metrics, and configuration together so
        # the only download boundary can be replaced without changing the repo.
        self.repo = self.root / 'local repo'
        shutil.copytree(coding_fixture.ROOT / 'src', self.repo / 'src')
        shutil.copytree(coding_fixture.ROOT / 'config', self.repo / 'config')
        for name in ('EVOSUBSTER_REUSE_FROM', 'EVOSUBSTER_REUSE_PLAN',
                     'EVOSUBSTER_REUSE_LAST_VERSION', 'LASTAL_THREADS',
                     'LASTAL_PIPELINE_MODE', 'EVOSUBSTER_LAST_OPTIONS',
                     'EVOSUBSTER_SELECTION_SOURCE', 'EVOSUBSTER_METADATA_INPUT'):
            self.env.pop(name, None)
        self.accessions = ['GCA_000000001.1', 'GCA_000000002.1', 'GCA_000000003.1']
        self.genomes = self.root / 'local assemblies'
        self.genomes.mkdir()
        rows = {}
        for accession, old, sequence in zip(self.accessions, self.fasta,
                                             ('ACGTNNNN', 'CCGTNNNN', 'GGGTNNNN')):
            name = Path(old).parent.name
            folder = self.genomes / name
            folder.mkdir()
            fasta = folder / 'genomic.fna'
            fasta.write_text('>chr1\n' + sequence + '\n')
            summary = folder / 'assembly_summary.json'
            summary.write_text(json.dumps({'reports': [{'accession': accession,
                                             'organism': {'organism_name': name}}]}))
            # Use the actual assembly-bundle contract with independently hashed
            # fixture bytes; no download/network operation creates these inputs.
            (folder / 'assembly_bundle.json').write_text(json.dumps({
                'schema_version': 2, 'accession': accession, 'fasta': 'genomic.fna',
                'fasta_sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(),
                'gff': None, 'gff_sha256': None}))
            rows[accession] = [name, str(fasta), str(summary), name, name, '']
        self.fasta = [rows[accession][1] for accession in self.accessions]
        self.download_rows = self.root / 'download rows.json'
        self.download_rows.write_text(json.dumps(rows))
        self.download_log = self.root / 'download calls.jsonl'
        self.training_log = self.root / 'training calls'
        self.scheduler_log = self.root / 'scheduler calls'
        self.env.update(LOCAL_DOWNLOAD_ROWS=str(self.download_rows),
                        LOCAL_DOWNLOAD_LOG=str(self.download_log),
                        LOCAL_TRAINING_LOG=str(self.training_log),
                        LOCAL_SCHEDULER_LOG=str(self.scheduler_log))
        (self.repo / 'src/dwl_organism.sh').write_text('''#!/bin/bash
python3 - "$1" <<'END'
import json, os, sys
with open(os.environ['LOCAL_DOWNLOAD_LOG'], 'a') as log:
    log.write(json.dumps({'accession': sys.argv[1],
                         'threads': os.environ.get('LASTAL_THREADS'),
                         'mode': os.environ.get('LASTAL_PIPELINE_MODE')}) + '\\n')
with open(os.environ['LOCAL_DOWNLOAD_ROWS']) as data:
    print('|'.join(json.load(data)[sys.argv[1]]))
END
''')
        training = self.align / 'last_train.sh'
        training.write_text(training.read_text().replace(
            '#!/bin/bash\n', '#!/bin/bash\nprintf "trained\\n" >> "$LOCAL_TRAINING_LOG"\n', 1))
        for name in ('qsub', 'qstat'):
            path = self.bin / name
            path.write_text('#!/bin/bash\nprintf "%s\\n" "$0" >> "$LOCAL_SCHEDULER_LOG"\nexit 99\n')
            path.chmod(0o755)
        self.trios = self.root / 'selected_trios.tsv'
        self.trios.write_text('out_acc\tin1_acc\tin2_acc\n' + '\t'.join(self.accessions) + '\n')

    def run_command(self, entrypoint, *args):
        return subprocess.run(['bash', str(self.repo / 'src' / entrypoint), *map(str, args)],
                              env=self.env, capture_output=True, text=True, timeout=20)

    def snapshot(self):
        return {str(path.relative_to(self.root)):
                (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
                if path.is_file() else None
                for path in self.root.rglob('*')}

    def test_trio_child_imports_native_stages_and_preview_has_no_side_effects(self):
        source_root = self.root / 'source results'
        first = self.run_command('sbst.sh', 'source-date', *self.fasta, 'NO_GFF_FILE',
                                 '--out-dir', source_root, '--thread', '2', '--idt-only',
                                 '--accession-ids', *self.accessions)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        source_manifest = next(source_root.rglob('metadata_manifest.json'))
        source_run = source_manifest.parent.parent
        original = json.loads(source_manifest.read_text())
        self.assertEqual(original['status'], 'identity_complete')
        expected = {'gc:org1', 'gc:org2', 'gc:org3',
                    'train:org1-org2', 'train:org1-org3', 'train:org2-org3'}
        self.assertEqual(set(original['stage_records']), expected)
        self.assertEqual(self.training_log.read_text().splitlines(), ['trained'] * 3)
        source_bytes = source_manifest.read_bytes()
        source_stamp = source_manifest.stat().st_mtime_ns

        target = self.root / 'new results'
        args = ['--trios', self.trios, 'new-date', '--genome-dir', self.genomes,
                '--out-dir', target, '--thread', '2', '--idt-only',
                '--reuse-from', source_run, '--lastal-threads', '1', '--alignment-mode', 'staged']
        before = self.snapshot()
        preview = self.run_command('sbst_fromDwl.sh', *args, '--reuse-plan')
        self.assertEqual(preview.returncode, 0, preview.stdout + preview.stderr)
        decisions = json.loads(preview.stdout)
        reusable = {row['kind'] + ':' + row['slot'] for row in decisions['stages']
                    if row['action'] == 'reuse'}
        self.assertTrue(expected.issubset(reusable), preview.stdout)
        self.assertEqual(self.snapshot(), before)
        self.assertFalse(target.exists())
        self.assertFalse(self.download_log.exists())
        self.assertFalse(self.scheduler_log.exists())

        result = self.run_command('sbst_fromDwl.sh', *args)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        manifests = list(target.rglob('metadata_manifest.json'))
        self.assertEqual(len(manifests), 1)
        child = json.loads(manifests[0].read_text())
        self.assertEqual(child['status'], 'identity_complete')
        self.assertEqual(child['date'], 'new-date')
        self.assertEqual(child['options']['last']['lastal'][:2], ['-P', '1'])
        self.assertEqual(child['options']['last']['last_train'][:2], ['-P', '2'])
        self.assertEqual(set(child['stage_records']), expected)
        for key in expected:
            with self.subTest(stage=key):
                record = child['stage_records'][key]
                self.assertEqual(record['provenance']['source_run'], str(source_run))
                self.assertEqual(record['provenance']['policy'], 'same')
                self.assertEqual(record['outputs'][0]['sha256'],
                                 original['stage_records'][key]['outputs'][0]['sha256'])
        self.assertEqual(self.training_log.read_text().splitlines(), ['trained'] * 3)
        calls = [json.loads(line) for line in self.download_log.read_text().splitlines()]
        self.assertEqual([row['accession'] for row in calls], self.accessions)
        self.assertTrue(all(row['threads'] == '1' and row['mode'] == 'staged' for row in calls))
        self.assertEqual(source_manifest.read_bytes(), source_bytes)
        self.assertEqual(source_manifest.stat().st_mtime_ns, source_stamp)
        self.assertFalse(self.scheduler_log.exists())

        # Previewing an existing destination must also leave imported outputs,
        # source files, logs, manifests, and directory structure untouched.
        before = self.snapshot()
        repeat = self.run_command('sbst_fromDwl.sh', *args, '--reuse-plan')
        self.assertEqual(repeat.returncode, 0, repeat.stdout + repeat.stderr)
        self.assertEqual(self.snapshot(), before)


if __name__ == '__main__':
    unittest.main()
