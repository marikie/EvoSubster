import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from test_collect_run_summary import create_run,script

ROOT=Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('Rscript'),'Rscript required')
class PCAManifestMainTests(unittest.TestCase):
    def test_main_uses_manifest_run_and_handles_constant_channels(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex',prefix='pca-main-v2-') as tmp:
            root=Path(tmp)
            run=create_run(root/'inputs'/'arbitrary-parent')
            manifest=script.runs.load(run)
            baseline=ROOT/'test/fixtures/pipeline_counts/20261031101732'
            for artifact in manifest['artifacts']:
                if artifact['kind']!='single': continue
                name='Oikalb2' if artifact['slot']=='org2' else 'Oikvan3'
                source=next(path for path in (baseline/name/'singlenuc').glob('*.tsv') if '_ncds' not in path.name)
                target=run/artifact['path']
                target.write_bytes(source.read_bytes())
                artifact['sha256']=script.runs.sha256(target)
            script.runs.atomic_json(run/script.runs.MANIFEST,manifest)
            result=subprocess.run(['Rscript',str(ROOT/'src/visualize/pca_2d.R'),str(root/'inputs'),
                '--tsv-pattern','*.tsv','--output-dir',str(root/'pca')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+'\n'+result.stderr)
            for variant in ('logratio','observed'):
                path=root/'pca'/f'clustering_results_{variant}.csv'
                self.assertTrue(path.is_file(),result.stderr)
                self.assertEqual(len(path.read_text().splitlines()),3)
