import csv
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from test_coding_region_pipeline import CodingRegionPipelineTest
ROOT = Path(__file__).resolve().parents[1]


class TrainingFailureTests(CodingRegionPipelineTest):
    def test_each_training_failure_reaches_idt_only_exit(self):
        original = (self.align / 'last_train.sh').read_text()
        for pair in ('Alpone1:Bettwo2', 'Alpone1:Gamthr3', 'Bettwo2:Gamthr3'):
            with self.subTest(pair=pair):
                self.out = self.root / ('failure-' + pair.replace(':', '-'))
                (self.align / 'last_train.sh').write_text(
                    '#!/bin/bash\nif [ "$4:$5" = "' + pair + '" ]; then exit 42; fi\n' + original)
                result = self.run_pipeline()
                self.assertNotEqual(result.returncode, 0)
                manifest = next(self.out.rglob('metadata_manifest.json'))
                self.assertNotIn('"status": "identity_complete"', manifest.read_text())

    def test_incomplete_successful_training_is_rejected(self):
        (self.align / 'last_train.sh').write_text('#!/bin/bash\nprintf "# substitution percent identity: 95\\n" > "$6/${4}2${5}_${1}.train"\n')
        result = self.run_pipeline()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('incomplete training output', result.stderr)


class PlotFailureTests(unittest.TestCase):
    def test_plot_failure_does_not_publish_figures_or_completion(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex') as tmp:
            work = Path(tmp)
            binary = work/'bin'
            binary.mkdir()
            stub=binary/'Rscript'
            stub.write_text('#!/bin/bash\nif [ "$1" = "-e" ]; then exit 0; fi\nexit 42\n')
            stub.chmod(0o755)
            inputs=[]
            for i in range(8):
                path=work/f'input{i}.tsv'
                path.write_text('fixture')
                inputs.append(str(path))
            for failed in ('-e','sbmut.R','logRatioPlot.R','dinucleotide-plot.R'):
                with self.subTest(failed=failed):
                    stub.write_text('#!/bin/bash\ncase "$1" in *"$FAIL_PLOT") exit 42 ;; esac\nexit 0\n')
                    result=subprocess.run(['bash',str(ROOT/'src/generate_graphs.sh'),*inputs,
                        str(ROOT/'src/visualize'),'figs/org2','figs/org3'],cwd=work,
                        env=dict(os.environ,FAIL_PLOT=failed,PATH=str(binary)+':'+os.environ['PATH']),capture_output=True)
                    self.assertNotEqual(result.returncode,0)
                    self.assertIn(b'non-zero exit status 42',result.stderr)
                    self.assertFalse(list(work.rglob('*.pdf')))
                    self.assertFalse((work/'metadata/plot_completion.json').exists())



@unittest.skipUnless(shutil.which('Rscript'),'Rscript required')
class DBSNullTests(unittest.TestCase):
    def run_counts(self, counts, duplicate=False, missing=False):
        rows=[('AT>CA',counts[0]),('AT>CC',counts[1]),('AT>CG',counts[2]),
              ('AT>GA',counts[3]),('AT>GC',counts[4]),('AT>TA',counts[5])]
        if duplicate: rows.append(rows[0])
        if missing: rows.pop()
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex') as tmp:
            path=Path(tmp)/'counts.tsv'
            with path.open('w') as handle:
                writer=csv.writer(handle,delimiter='\t')
                writer.writerow(['mutType','mutNum','totalRootNum'])
                writer.writerows((channel,count,10000) for channel,count in rows)
            return subprocess.run(['Rscript',str(ROOT/'src/statistics/chi2_context_dinuc.R'),str(path)],text=True,capture_output=True)

    def test_low_counts_are_reproducible(self):
        first=self.run_counts([2,0,1,0,0,0]); second=self.run_counts([2,0,1,0,0,0])
        self.assertEqual(first.returncode,0,first.stderr)
        self.assertEqual(first.stdout,second.stdout)
        self.assertIn('multiMC',first.stdout)

    def test_zero_counts_skip_p_value(self):
        result=self.run_counts([0]*6)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('skip',result.stdout)

    def test_missing_or_duplicate_channel_is_rejected(self):
        for options in ({'missing':True},{'duplicate':True}):
            result=self.run_counts([1]*6,**options)
            self.assertNotEqual(result.returncode,0)

    def test_uniform_original_nine_are_uniform_after_folding(self):
        # AT has three paired and three singleton target categories.
        rows=[('AT>CA',200),('AT>CC',200),('AT>CG',100),
              ('AT>GA',200),('AT>GC',100),('AT>TA',100)]
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex') as tmp:
            path=Path(tmp)/'counts.tsv'
            with path.open('w') as handle:
                writer=csv.writer(handle,delimiter='\t')
                writer.writerow(['mutType','mutNum','totalRootNum'])
                writer.writerows((channel,count,10000) for channel,count in rows)
            result=subprocess.run(['Rscript',str(ROOT/'src/statistics/chi2_context_dinuc.R'),str(path)],text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            line=next(line for line in result.stdout.splitlines() if line.startswith('AT '))
            self.assertEqual(float(line.split()[4]),1.0)
            self.assertIn('B=10000, seed=1',result.stdout)
