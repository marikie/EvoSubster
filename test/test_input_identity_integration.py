import json
from pathlib import Path
import unittest
from test_coding_region_pipeline import CodingRegionPipelineTest


class InputIdentityIntegrationTests(CodingRegionPipelineTest):
    def test_same_labels_date_use_content_order_and_options(self):
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
            result=self.run_pipeline(*options)
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
