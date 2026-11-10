"""CLI contracts: previews never publish outputs; successful stages survive interruption."""
import json
import shutil
from pathlib import Path
import subprocess
from test_coding_region_pipeline import CodingRegionPipelineTest, ROOT

class ReuseIntegrationTests(CodingRegionPipelineTest):
    def test_preview_does_not_create_output_or_execute_training(self):
        result=self.run_pipeline('--reuse-plan')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('recompute',result.stdout)
        self.assertFalse(self.out.exists())

    def test_unavailable_download_input_preview_never_downloads(self):
        genomes=self.root/'not-downloaded'
        result=subprocess.run(['bash',str(ROOT/'src/sbst_fromDwl.sh'),'20261110',
            'GCA_000000001.1','GCA_000000002.1','GCA_000000003.1',
            '--genome-dir',str(genomes),'--out-dir',str(self.out),'--reuse-plan'],
            env=self.env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('unavailable',result.stdout)
        self.assertFalse(genomes.exists())
        self.assertFalse(self.out.exists())

    def test_successful_metrics_registered_before_training_failure(self):
        (self.align/'last_train.sh').write_text('#!/bin/bash\nexit 17\n')
        result=self.run_pipeline()
        self.assertNotEqual(result.returncode,0)
        value=json.loads(next(self.out.rglob('metadata_manifest.json')).read_text())
        self.assertIn('gc:org1',value.get('stage_records',{}))
        gc=next(self.misc.glob('*gcContent*'))
        stamp=gc.stat().st_mtime_ns
        result=self.run_pipeline()
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(gc.stat().st_mtime_ns,stamp)

    def test_cross_run_reuse_preserves_files_under_new_date_names(self):
        result=self.run_pipeline()
        self.assertEqual(result.returncode,0,result.stderr)
        source=next(self.out.rglob('metadata_manifest.json')).parent.parent
        source_manifest=(source/'metadata/metadata_manifest.json').read_bytes()
        target=self.root/'different-results'
        command=['bash',str(ROOT/'src/sbst.sh'),'20261110',*self.fasta,'NO_GFF_FILE',
            '--out-dir',str(target),'--idt-only','--reuse-from',str(source)]
        result=subprocess.run(command,env=self.env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        value=json.loads(next(target.rglob('metadata_manifest.json')).read_text())
        self.assertTrue(value['stage_records']['gc:org1'].get('provenance'))
        self.assertEqual(source_manifest,(source/'metadata/metadata_manifest.json').read_bytes())
        self.assertTrue(list(target.rglob('*gcContent_20261110.out')))

    def test_changed_gc_implementation_recomputes_without_retraining(self):
        metrics=self.root/'metrics'
        metrics.mkdir()
        shutil.copy2(ROOT/'src/metrics/gc_content.sh',metrics/'gc_content.sh')
        self.env['METRICS_DIR_OVERRIDE']=str(metrics)
        first=self.run_pipeline()
        self.assertEqual(first.returncode,0,first.stderr)
        run=next(self.out.rglob('metadata_manifest.json')).parent.parent
        trains={str(p):p.stat().st_mtime_ns for p in run.rglob('*.train')}
        (metrics/'gc_content.sh').write_text("#!/bin/bash\nprintf 'Total GC content: 64.00%%\\n'\n")
        second=self.run_pipeline()
        self.assertEqual(second.returncode,0,second.stderr)
        self.assertEqual(next(self.misc.glob('*gcContent*')).read_text(),'Total GC content: 64.00%\n')
        self.assertEqual(trains,{str(p):p.stat().st_mtime_ns for p in run.rglob('*.train')})

    def test_preview_never_imports_implicitly_from_other_runs(self):
        first=self.run_pipeline()
        self.assertEqual(first.returncode,0,first.stderr)
        target=next(self.out.rglob('metadata_manifest.json')).parent.parent
        other=subprocess.run(['bash',str(ROOT/'src/sbst.sh'),'20261110',*self.fasta,
            'NO_GFF_FILE','--out-dir',str(self.out),'--idt-only'],env=self.env,capture_output=True,text=True)
        self.assertEqual(other.returncode,0,other.stderr)
        target_manifest=json.loads((target/'metadata/metadata_manifest.json').read_text())
        outgroup=next(o['short_name'] for o in target_manifest['organisms'] if o['slot']=='org1')
        missing=next((target/'statistics/misc').glob(f'{outgroup}_gcContent*'))
        missing.unlink()
        before=(target/'metadata/metadata_manifest.json').read_bytes()
        preview=self.run_pipeline('--reuse-plan')
        self.assertEqual(preview.returncode,0,preview.stderr)
        value=json.loads(preview.stdout)
        row=next(r for r in value['stages'] if r['kind']=='gc' and r['slot']=='org1')
        self.assertEqual(row['action'],'recompute')
        self.assertEqual(before,(target/'metadata/metadata_manifest.json').read_bytes())
        self.assertFalse(missing.exists())

    def test_tree_preview_does_not_use_stale_selection_outputs(self):
        saved=self.out/'trio_selection'
        saved.mkdir(parents=True)
        (saved/'selected_trios.tsv').write_text('obsolete selection output')
        tree=self.root/'changed-tree.nwk'
        tree.write_text('(New_species,Other_species);\n')
        before=(saved/'selected_trios.tsv').read_bytes()
        result=subprocess.run(['bash',str(ROOT/'src/sbst_fromDwl.sh'),'--tree',str(tree),
            '20261110','--out-dir',str(self.out),'--reuse-plan'],env=self.env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('Tree selection',result.stdout)
        self.assertEqual(before,(saved/'selected_trios.tsv').read_bytes())
