import csv
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
            self.assertFalse(list((root/'pca').glob('*_filtered*')))
            settings=json.loads((root/'pca/pca_run_settings.json').read_text())
            self.assertEqual(settings['status'],'completed')
            self.assertEqual(settings['identity_source'],'last-train-penultimate')
            self.assertEqual(settings['usable_samples'],2)
            self.assertIsNone(settings['effective_thresholds']['ingroup_idt_threshold'])
            with (root/'pca/pca_trio_filter.csv').open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['usable_samples'],'2')
            self.assertEqual(rows[0]['identity_source'],'last-train-penultimate')
            for variant in ('logratio','observed'):
                path=root/'pca'/f'clustering_results_{variant}.csv'
                self.assertTrue(path.is_file(),result.stderr)
                self.assertEqual(len(path.read_text().splitlines()),3)

    def test_all_excluded_records_saved(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex',prefix='pca-excluded-') as tmp:
            root=Path(tmp)
            create_run(root/'inputs'/'parent')
            result=subprocess.run(['Rscript',str(ROOT/'src/visualize/pca_2d.R'),str(root/'inputs'),
                '--ingroup-idt-threshold','98','--output-dir',str(root/'pca')],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            settings=json.loads((root/'pca/pca_run_settings.json').read_text())
            self.assertEqual(settings['status'],'failed')
            self.assertEqual(settings['usable_samples'],0)
            with (root/'pca/pca_trio_filter.csv').open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(rows[0]['included'],'FALSE')
            self.assertIn('failed',rows[0]['reason'])

    def test_active_filter_selects_trio_and_cli_override_order(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex',prefix='pca-selection-') as tmp:
            root=Path(tmp)
            baseline=ROOT/'test/fixtures/pipeline_counts/20261031101732'
            for key in ('retained','excluded'):
                run=create_run(root/'inputs'/key)
                manifest=script.runs.load(run)
                manifest['trio_key']=key
                for artifact in manifest['artifacts']:
                    if artifact['kind']=='single':
                        name='Oikalb2' if artifact['slot']=='org2' else 'Oikvan3'
                        source=next(p for p in (baseline/name/'singlenuc').glob('*.tsv') if '_ncds' not in p.name)
                        target=run/artifact['path']
                        target.write_bytes(source.read_bytes())
                        artifact['sha256']=script.runs.sha256(target)
                if key=='excluded':
                    manifest['organisms'][2]['raw_organism_name']='Different gamma'
                script.runs.atomic_json(run/script.runs.MANIFEST,manifest)
            for index,options in enumerate([
                ['--idt-threshold','99','--ingroup-idt-threshold','96','--outgroup-idt-threshold','89'],
                ['--outgroup-idt-threshold','89','--ingroup-idt-threshold','96','--idt-threshold','99']]):
                output=root/f'pca{index}'
                result=subprocess.run(['Rscript',str(ROOT/'src/visualize/pca_2d.R'),str(root/'inputs'),
                    '--tsv-pattern','*.tsv','--output-dir',str(output),
                    '--min-idt-gap-between-out-in','10',*options],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                settings=json.loads((output/'pca_run_settings.json').read_text())
                self.assertEqual(settings['included_trios'],1)
                self.assertEqual(settings['usable_samples'],2)
                self.assertEqual(settings['effective_thresholds']['ingroup_idt_threshold'],96)
                self.assertEqual(settings['effective_thresholds']['outgroup_idt_threshold'],89)
                with (output/'pca_trio_filter.csv').open() as handle:
                    rows={row['trio_id']:row for row in csv.DictReader(handle)}
                self.assertEqual(rows['retained']['min_idt_gap_between_out_in'],'waived')
                self.assertEqual(rows['retained']['usable_samples'],'2')
                self.assertEqual(rows['excluded']['included'],'FALSE')
                self.assertEqual(rows['excluded']['min_idt_gap_between_out_in'],'failed')
                for variant in ('logratio','observed'):
                    text=(output/f'clustering_results_{variant}.csv').read_text()
                    self.assertIn('retained',text)
                    self.assertNotIn('excluded',text)

    def test_legacy_inputs_filter_audit_and_missing_sample_reason(self):
        from test_pca_input_adapter import create_legacy
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex',prefix='pca-legacy-main-') as tmp:
            root = Path(tmp)
            retained = create_legacy(root/'inputs')
            create_legacy(root/'inputs', name='Alt1_Alt2_Alt3', nested=False)
            # One selected sample missing: must remain visible in trio audit.
            missing = next((retained/'statistics/In3/singlenuc').glob('*_20251027.tsv'))
            missing.unlink()
            command = ['Rscript',str(ROOT/'src/visualize/pca_2d.R'),str(root/'inputs'),
                       '--tsv-pattern','*.tsv','--output-dir',str(root/'pca'),
                       '--idt-threshold','80','--ingroup-idt-threshold','80',
                       '--outgroup-idt-threshold','80','--min-idt-gap-between-out-in','0']
            result = subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            settings = json.loads((root/'pca/pca_run_settings.json').read_text())
            self.assertEqual(settings['trios'],2)
            self.assertEqual(settings['included_trios'],2)
            self.assertEqual(settings['usable_samples'],3)
            with (root/'pca/pca_trio_filter.csv').open() as handle:
                rows = {row['trio_id']:row for row in csv.DictReader(handle)}
            self.assertEqual(rows['Out1_In2_In3']['usable_samples'],'1')
            self.assertIn('org3: Missing canonical single TSV',rows['Out1_In2_In3']['reason'])
            self.assertEqual(rows['Alt1_Alt2_Alt3']['usable_samples'],'2')
            self.assertFalse(list((root/'pca').glob('*_filtered*')))
