import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('jq'),'jq required')
class AssemblyBundleDownloadTests(unittest.TestCase):
    def test_assembly_isolation_resolution_and_ownership(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp/codex',prefix='bundle-download-') as tmp:
            work=Path(tmp)
            fixtures=work/'fixtures'; fixtures.mkdir()
            for accession,sequence,annotation in [('GCA_000000001.1','ACGT',True),('GCA_000000002.1','AAAA',False),('GCA_000000001.2','AGGT',True)]:
                with zipfile.ZipFile(fixtures/(accession+'.zip'),'w') as archive:
                    prefix='ncbi_dataset/data/'+accession+'/'
                    archive.writestr(prefix+accession+'_genomic.fna','>chr1\n'+sequence+'\n')
                    if annotation: archive.writestr(prefix+'genomic.gff','chr1\ttest\tCDS\t1\t3\t.\t+\t0\tID='+accession+'\n')
            binary=work/'bin'; binary.mkdir()
            curl=binary/'curl'
            curl.write_text('''#!/usr/bin/env python3
import json,os,sys,shutil
from pathlib import Path
args=sys.argv[1:]
url=next(x for x in args if x.startswith('http'))
requested=url.split('/accession/')[1].split('/')[0]
accession=requested if '.' in requested else requested+'.'+os.environ.get('LATEST_VERSION','1')
if '-o' in args: shutil.copyfile(Path(os.environ['FIXTURE_ROOT'])/(accession+'.zip'),args[args.index('-o')+1])
else: print(json.dumps({'reports':[{'accession':accession,'organism':{'organism_name':'Test species','tax_id':123}}]}))
''')
            curl.chmod(0o755)
            genomes=work/'genomes'; genomes.mkdir()
            legacy=genomes/'Test_species'; legacy.mkdir(); (legacy/'valuable').write_text('preserve')
            env=dict(os.environ,PATH=str(binary)+':'+os.environ['PATH'],FIXTURE_ROOT=str(fixtures))
            def download(accession,extra=(),latest='1'):
                result=subprocess.run(['bash',str(ROOT/'src/dwl_organism.sh'),accession,'--out-dir',str(genomes),'--no-taxonomy',*extra],
                    env=dict(env,LATEST_VERSION=latest),capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                fields=result.stdout.strip().splitlines()[-1].split('|')
                return Path(fields[1])
            ownership=work/'owned.tsv'
            first=download('GCA_000000001.1',('--artifact-manifest',str(ownership)))
            rows=list(csv.DictReader(ownership.read_text().splitlines(),delimiter='\t'))
            self.assertEqual(rows[0]['relative_path'],'assemblies/GCA_000000001.1')
            self.assertTrue((genomes/rows[0]['relative_path']).is_dir())
            second=download('GCA_000000002.1')
            self.assertNotEqual(first,second)
            self.assertIsNone(json.loads((second.parent/'assembly_bundle.json').read_text())['gff'])
            again=download('GCA_000000001.1')
            self.assertEqual(first,again)
            self.assertIn('GCA_000000001.1',(again.parent/'genomic.gff').read_text())
            self.assertEqual(download('GCA_000000001'),first)
            newer=download('GCA_000000001',latest='2')
            self.assertIn('GCA_000000001.2',str(newer))
            self.assertNotEqual(first,newer)
            existing_owned=work/'existing-owned.tsv'
            download('GCA_000000001.1',('--force','--artifact-manifest',str(existing_owned)))
            existing=list(csv.DictReader(existing_owned.read_text().splitlines(),delimiter='\t'))
            self.assertTrue(existing)
            for item in existing:
                self.assertTrue(item['relative_path'].startswith('assemblies/GCA_000000001.1/'))
            # New-accession cleanup targets the owned accession bundle, never the legacy organism folder.
            shutil.rmtree(genomes/rows[0]['relative_path'])
            self.assertEqual((legacy/'valuable').read_text(),'preserve')
            self.assertTrue(second.is_file())
