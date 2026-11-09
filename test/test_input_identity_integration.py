import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time
import unittest
from test_coding_region_pipeline import CodingRegionPipelineTest, ROOT


class InputIdentityIntegrationTests(CodingRegionPipelineTest):
    def test_readable_metadata_and_renamed_old_run_resume(self):
        metadata = self.root / 'metadata.json'
        entries = [dict(raw_organism_name=name, short_name=short)
                   for name, short in zip(('Alpha full species', 'Beta full species', 'Gamma full species'),
                                          ('Alpone1', 'Bettwo2', 'Gamthr3'))]
        metadata.write_text(json.dumps(entries))
        self.env['EVOSUBSTER_METADATA_INPUT'] = str(metadata)
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stderr)
        path = next(self.out.glob('*/*/metadata/metadata_manifest.json'))
        value = json.loads(path.read_text())
        run = path.parent.parent
        self.assertIn('Alpha_full_species_local-', run.parent.name)
        old = self.out / value['trio_key'] / ('20260101-' + value['run_key'])
        old.parent.mkdir()
        run.rename(old)
        run.parent.rmdir()
        for index, fasta in enumerate(self.fasta):
            moved = self.root / f'renamed_{index}'
            Path(fasta).parent.rename(moved)
            self.fasta[index] = str(moved / 'genomic.fna')
            entries[index]['raw_organism_name'] = 'Changed display name'
        metadata.write_text(json.dumps(entries))
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(self.out.glob('*/*/metadata/metadata_manifest.json')), [old / 'metadata/metadata_manifest.json'])
        self.assertEqual(json.loads((old / 'metadata/metadata_manifest.json').read_text())['created_at'], value['created_at'])

    def test_allocation_serialization_collision_and_active_run_lock(self):
        # Force a cryptographic prefix collision at the descriptor boundary;
        # execute the real production allocation/prepare/flock path in sbst.sh.
        wrapper = self.bin / 'python3'
        wrapper.write_text(f'''#!{sys.executable}
import os, sys, time
if len(sys.argv) > 1 and sys.argv[1].endswith('/run_manifest.py') and os.environ.get('FORCED_RUN_KEY'):
    sys.path.insert(0, {str(ROOT / 'src/lib')!r})
    import run_manifest
    original = run_manifest.descriptor
    def descriptor(args):
        value = original(args)
        value['run_key'] = os.environ['FORCED_RUN_KEY']
        return value
    run_manifest.descriptor = descriptor
    if len(sys.argv) > 2 and sys.argv[2] == 'prepare' and os.environ.get('HOLD_PREPARE'):
        open(os.environ['PREPARE_STARTED'], 'w').close()
        while not os.path.exists(os.environ['PREPARE_RELEASE']):
            time.sleep(0.05)
    sys.argv = sys.argv[1:]
    run_manifest.main()
else:
    os.execv({sys.executable!r}, [{sys.executable!r}, *sys.argv[1:]])
''')
        wrapper.chmod(0o755)
        started, release = self.root / 'started', self.root / 'release'
        preparing, publish = self.root / 'preparing', self.root / 'publish'
        profile = self.align / 'last_train.sh'
        profile.write_text('''#!/bin/bash
if [ -n "${HOLD_TRAINING:-}" ]; then
    touch "$TRAIN_STARTED"
    while [ ! -e "$TRAIN_RELEASE" ]; do sleep 0.05; done
fi
''' + profile.read_text())
        command = ['bash', str(ROOT / 'src/sbst.sh'), '20260101', *self.fasta,
                   'NO_GFF_FILE', '--out-dir', str(self.out), '--idt-only']
        first_env = dict(self.env, FORCED_RUN_KEY='a' * 64, HOLD_TRAINING='1',
                         TRAIN_STARTED=str(started), TRAIN_RELEASE=str(release),
                         HOLD_PREPARE='1', PREPARE_STARTED=str(preparing), PREPARE_RELEASE=str(publish))
        first = subprocess.Popen(command, env=first_env, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True)
        second = None
        try:
            deadline = time.monotonic() + 30
            while not preparing.exists() and first.poll() is None and time.monotonic() < deadline:
                time.sleep(0.05)
            self.assertTrue(preparing.exists(), 'First run never reached prepare')
            with (self.out / '.run-allocation.lock').open('a') as allocation:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(allocation, fcntl.LOCK_EX | fcntl.LOCK_NB)
            second = subprocess.Popen(command, env=dict(self.env, FORCED_RUN_KEY='a' * 12 + 'b' * 52),
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # During the unpublished allocation window the second process must wait.
            time.sleep(0.3)
            self.assertIsNone(second.poll())
            self.assertFalse(list(self.out.glob('*/*/metadata/metadata_manifest.json')))
            publish.touch()
            deadline = time.monotonic() + 30
            while not started.exists() and first.poll() is None and time.monotonic() < deadline:
                time.sleep(0.05)
            self.assertTrue(started.exists(), 'First run never reached training')
            same = subprocess.run(command, env=dict(self.env, FORCED_RUN_KEY='a' * 64),
                                  capture_output=True, text=True, timeout=30)
            self.assertNotEqual(same.returncode, 0)
            self.assertIn('run is already active', same.stderr)
            other_stdout, other_stderr = second.communicate(timeout=30)
            self.assertEqual(second.returncode, 0, other_stderr)
            runs = list(self.out.glob('*/*/metadata/metadata_manifest.json'))
            self.assertEqual({p.parent.parent.name for p in runs},
                             {'20260101-' + 'a' * 12, '20260101-' + 'a' * 12 + 'b' * 4})
        finally:
            publish.touch()
            release.touch()
            stdout, stderr = first.communicate(timeout=30)
            if second is not None and second.poll() is None:
                second.communicate(timeout=30)
        self.assertEqual(first.returncode, 0, stderr)

    def test_same_labels_date_use_content_order_and_options(self):
        accession_options = ('--accession-ids', 'GCA_000000001.1', 'GCA_000000002.1', 'GCA_000000003.1')
        profile=self.align/'last_train.sh'
        original=profile.read_text()
        log=self.root/'training-calls.txt'
        self.env['TRAIN_CALL_LOG']=str(log)
        profile.write_text('''#!/bin/bash
output="$6/${4}2${5}_${1}.train"
if [ -s "$output" ]; then exit 0; fi
printf 'trained\\n' >> "$TRAIN_CALL_LOG"
'''+original)
        def run(*options):
            result=self.run_pipeline(*accession_options, *options)
            self.assertEqual(result.returncode,0,result.stderr)
            return {json.loads(p.read_text())['run_key'] for p in self.out.rglob('metadata_manifest.json')}
        self.assertEqual(len(run()),1)
        self.assertEqual(len(log.read_text().splitlines()),3)
        self.assertEqual(len(run()),1)
        self.assertEqual(len(log.read_text().splitlines()),3)
        Path(self.fasta[0]).write_text('>chr1\nAAAANNNN\n')
        self.assertEqual(len(run()),2)
        self.assertEqual(len(log.read_text().splitlines()),6)
        self.fasta[0],self.fasta[1]=self.fasta[1],self.fasta[0]
        self.assertEqual(len(run()),3)
        self.assertEqual(len(run('--lastal-options=-H2')),4)
        self.annotate(0)
        result = self.run_pipeline(*accession_options, outgroup=self.gffs[0])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len({json.loads(p.read_text())['run_key'] for p in self.out.rglob('metadata_manifest.json')}), 5)
        self.gffs[0].write_text('chr1\ttest\tCDS\t1\t4\t.\t+\t0\tID=cds\n')
        result = self.run_pipeline(*accession_options, outgroup=self.gffs[0])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len({json.loads(p.read_text())['run_key'] for p in self.out.rglob('metadata_manifest.json')}), 6)
