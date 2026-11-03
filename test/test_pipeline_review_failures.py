import csv
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


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
