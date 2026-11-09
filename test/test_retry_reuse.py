import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
import run_manifest as manifest
import retry_reuse as reuse
import last_options

TRAIN = '''# substitution percent identity: 95.0
#last -t 1
#last -a 1
#last -A 1
#last -b 1
#last -B 1
#last -S 1
# score matrix (query letters = columns, reference letters = rows):
A 1 -1 -1 -1
C -1 1 -1 -1
G -1 -1 1 -1
T -1 -1 -1 1
'''
MAF = '##maf version=1\n\na score=4\ns ref 0 4 + 10 ACGT\ns query 1 4 - 12 ACGT\n\n'
TIMING = '\nreal\t0m1.234s\nuser\t0m1.000s\nsys\t0m0.234s\n'


class ReuseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source, self.target = self.root/'source', self.root/'target'
        last = dict(schema=1, lastdb_train=['-P','48','-c','-u','RY4'],
                    lastdb_align=['-P','48','-c','-u','RY128'],
                    last_train=['-P','48','--revsym','-C','2'],
                    lastal=['-P','48','-i','1M','-H','1','-C','2','--split-f','MAF+','--split-b','8T'],
                    last_split=['-r'])
        accessions = ['GCF_040869185.1','GCA_026956025.1','GCF_049350105.2']
        organisms = []
        for i, accession in enumerate(accessions, 1):
            fasta = self.root / (accession + '.fa')
            fasta.write_text('>seq\nACGT\n')
            organisms.append(dict(slot=f'org{i}',accession=accession,short_name=f'Sp{i}',
                                  fasta_path=str(fasta), fasta_sha256=manifest.sha256(fasta),
                                  gff_path=None,gff_sha256=None))
        value = dict(schema_version=2,run_key='old',trio_key='trio',date='20261108',created_at='now',
                     options=dict(last=last,threads=48,last_version='lastal 1648',algorithm='v2',config={'sbst':'abc'}),
                     organisms=organisms,artifacts=[],status='pending')
        manifest.prepare(self.source, copy.deepcopy(value))
        self.maf = 'intermediateFiles/Sp12Sp2_one2one_20261108.maf'
        self.train = 'intermediateFiles/Sp12Sp2_20261108.train'
        (self.source/'intermediateFiles').mkdir()
        (self.source/self.train).write_text(TRAIN)
        manifest.seal_file(self.source/self.train)
        (self.source/'statistics').mkdir()
        (self.source/'statistics/gc.out').write_text('0.5\n')
        manifest.checkpoint(self.source,'identity_complete',[
            f'train:org1-org2:{self.train}', 'gc:org1:statistics/gc.out'])
        # Legacy pair alignment was created after the identity checkpoint.
        (self.source/self.maf).write_text(MAF)
        (self.source/(self.maf+'.raw')).write_text('partial raw')
        (self.source/'intermediateFiles/old-db').mkdir()
        (self.source/'intermediateFiles/old-db/db.prj').write_text('legacy database')
        value['run_key']='new'
        value['options']['last_version']='lastal 1654'
        value['options']['last']['lastal'][1]='1'
        manifest.prepare(self.target,value)
        db='intermediateFiles/Sp1db_20261108.seed-RY128/Sp1db_20261108.seed-RY128'
        self.command = (f'time lastal -p {self.train} -P48 -i 1M -H 1 -C 2 --split-f=MAF+ --split-b=8T '
                        f'{db} {organisms[1]["fasta_path"]} | last-split -r >{self.maf}.raw')
        self.log=self.root/'job.log'
        self.log.write_text(f'pwd: {self.source}\n{self.command}\n'+TIMING+
                            f'time maf-linked {self.maf}.raw >{self.maf}\n'+TIMING)

    def test_import_copies_only_certified_files_and_preserves_source(self):
        before={str(p.relative_to(self.source)):manifest.sha256(p) for p in self.source.rglob('*') if p.is_file()}
        reuse.import_run(self.target,self.source,self.log)
        self.assertTrue(reuse.pair_ready(self.target,'org1-org2'))
        self.assertFalse(reuse.pair_ready(self.target,'org1-org3'))
        self.assertEqual((self.target/self.maf).read_text(),MAF)
        self.assertNotEqual((self.target/self.maf).stat().st_ino,(self.source/self.maf).stat().st_ino)
        self.assertFalse((self.target/(self.maf+'.raw')).exists())
        self.assertFalse((self.target/'intermediateFiles/old-db').exists())
        self.assertEqual(before,{str(p.relative_to(self.source)):manifest.sha256(p) for p in self.source.rglob('*') if p.is_file()})
        record=manifest.load(self.target)['legacy_reuse']
        self.assertEqual(record['source_options']['last_version'],'lastal 1648')
        manifest.prepare(self.target,manifest.load(self.target))
        self.assertTrue(reuse.pair_ready(self.target,'org1-org2'))
        before=manifest.sha256(self.target/manifest.MANIFEST)
        reuse.import_run(self.target,self.source,self.log)
        self.assertEqual(before,manifest.sha256(self.target/manifest.MANIFEST))

    def test_input_mismatch_refuses_before_copy(self):
        value=manifest.load(self.target)
        value['organisms'][0]['fasta_sha256']='different'
        manifest.atomic_json(self.target/manifest.MANIFEST,value)
        with self.assertRaisesRegex(ValueError,'inputs'):
            reuse.import_run(self.target,self.source,self.log)
        self.assertFalse((self.target/'statistics/gc.out').exists())

    def test_partial_maf_not_certified(self):
        (self.source/self.maf).write_text(MAF+'a score=5\ns ref 0 1 + 10 A\n')
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse(reuse.pair_ready(self.target,'org1-org2'))
        self.assertFalse((self.target/self.maf).exists())
        self.assertTrue((self.target/'statistics/gc.out').exists())

    def test_failed_or_incomplete_log_not_certified(self):
        for text in (self.log.read_text().replace(TIMING,''),
                     self.log.read_text().replace(self.command,self.command+'\nKilled'),
                     self.log.read_text().replace('-P48 -i 1M','-P48 -i 2M')):
            with self.subTest(text=text):
                self.log.write_text(text)
                reuse.import_run(self.target,self.source,self.log)
                self.assertFalse(reuse.pair_ready(self.target,'org1-org2'))
                value=manifest.load(self.target); value.pop('legacy_reuse',None)
                manifest.atomic_json(self.target/manifest.MANIFEST,value)

    def test_bad_training_or_stats_hash_is_never_blessed(self):
        (self.source/self.train).write_text(TRAIN+'tampered\n')
        (self.source/'statistics/gc.out').write_text('changed')
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse((self.target/self.maf).exists())
        self.assertFalse((self.target/'statistics/gc.out').exists())

    def test_policy_only_permits_declared_version_transition(self):
        for version in ('lastal 1655','lastal 1648'):
            value=manifest.load(self.target); value['options']['last_version']=version
            manifest.atomic_json(self.target/manifest.MANIFEST,value)
            with self.assertRaisesRegex(ValueError,'1648.*1654'):
                reuse.import_run(self.target,self.source,self.log)

    def test_tampered_import_fails_readiness(self):
        reuse.import_run(self.target,self.source,self.log)
        (self.target/self.maf).write_text(MAF+'# altered\n')
        self.assertFalse(reuse.pair_ready(self.target,'org1-org2'))

    def test_sorted_legacy_requires_completed_sort_and_refresh_preserves_provenance(self):
        reuse.import_run(self.target,self.source,self.log)
        sorted_name=self.maf[:-4]+'_sorted.maf'
        (self.target/self.maf).rename(self.target/sorted_name)
        reuse.refresh_sorted(self.target)
        self.assertTrue(reuse.pair_ready(self.target,'org1-org2'))
        record=manifest.load(self.target)['legacy_reuse']['pairs']['org1-org2']
        self.assertEqual(record['source_maf'],self.maf)
        self.assertEqual(record['maf'],sorted_name)
        manifest.prepare(self.target,manifest.load(self.target))
        reuse.import_run(self.target,self.source,self.log)
        self.assertTrue(reuse.pair_ready(self.target,'org1-org2'))

    def test_native_sorted_pair_requires_checkpoint_and_settings(self):
        sorted_name=self.maf[:-4]+'_sorted.maf'
        (self.target/'intermediateFiles').mkdir()
        (self.target/sorted_name).write_text(MAF)
        (self.target/self.train).write_text(TRAIN)
        manifest.seal_file(self.target/self.train)
        options_path=self.target/(self.maf+'.last_options.json')
        options_path.write_text(json.dumps(manifest.load(self.target)['options']['last']))
        self.assertFalse(reuse.current_pair_ready(self.target,'org1-org2'))
        manifest.checkpoint(self.target,'identity_complete',[f'train:org1-org2:{self.train}'])
        self.assertTrue(reuse.current_pair_ready(self.target,'org1-org2'))
        options_path.write_text('{}')
        self.assertFalse(reuse.current_pair_ready(self.target,'org1-org2'))

    def test_semantic_options_mismatch_refused(self):
        value=manifest.load(self.target)
        value['options']['last']['lastal'][3]='2M'
        manifest.atomic_json(self.target/manifest.MANIFEST,value)
        with self.assertRaisesRegex(ValueError,'semantic'):
            reuse.import_run(self.target,self.source,self.log)

    def test_alignment_reuse_is_limited_to_explicit_trios(self):
        for run in (self.source,self.target):
            value=manifest.load(run)
            value['organisms'][0]['accession']='GCF_000000001.1'
            manifest.atomic_json(run/manifest.MANIFEST,value)
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse((self.target/self.maf).exists())
        self.assertTrue((self.target/'statistics/gc.out').exists())

    def test_prepare_recovers_only_sealed_matching_stage_and_keeps_diagnostics(self):
        (self.target/'intermediateFiles').mkdir()
        (self.target/self.train).write_text(TRAIN)
        manifest.seal_file(self.target/self.train)
        manifest.checkpoint(self.target,'identity_complete',[f'train:org1-org2:{self.train}'])
        stage=self.target/(self.maf+'.lastal.maf')
        stage.write_text(MAF)
        diagnostic=self.target/(self.maf+'.last-split.stderr.123')
        diagnostic.write_text('Killed by scheduler')
        partial=self.target/(self.maf+'.raw.tmp.123');partial.write_text('partial')
        env=dict(THREAD_NUM='48',LASTAL_THREADS='1',LASTAL_PIPELINE_MODE='staged',EVOSUBSTER_LAST_OPTIONS='{}')
        with patch.dict(os.environ,env), patch('subprocess.check_output',return_value='lastal 1654'):
            value=manifest.load(self.target)
            inputs=[o['fasta_path'] for o in value['organisms'][:2]]+[str(self.target/self.train)]
            last_options.output_seal(stage,inputs)
            manifest.prepare(self.target,value)
            self.assertTrue(stage.exists())
            self.assertFalse(partial.exists())
            archived=self.target/'metadata/diagnostics'/diagnostic.name
            self.assertEqual(archived.read_text(),'Killed by scheduler')
            stage.write_text('tampered')
            manifest.prepare(self.target,manifest.load(self.target))
            self.assertFalse(stage.exists())

    def test_checkpoint_drops_consumed_stage_artifacts(self):
        (self.target/'intermediateFiles').mkdir()
        stage=self.target/(self.maf+'.lastal.maf');stage.write_text(MAF)
        manifest.checkpoint(self.target,'identity_complete',[])
        stage.unlink()
        (self.target/self.maf).write_text(MAF)
        manifest.checkpoint(self.target,'identity_complete',[])
        value=manifest.validated(self.target)
        self.assertNotIn(self.maf+'.lastal.maf',[a['path'] for a in value['artifacts']])

    def test_source_sorted_requires_own_successful_sort(self):
        sorted_name=self.maf[:-4]+'_sorted.maf'
        (self.source/self.maf).rename(self.source/sorted_name)
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse(reuse.pair_ready(self.target,'org1-org2'))
        value=manifest.load(self.target);value.pop('legacy_reuse')
        manifest.atomic_json(self.target/manifest.MANIFEST,value)
        self.log.write_text(self.log.read_text()+f'time maf-sort {self.maf} >{sorted_name}\n'+TIMING)
        reuse.import_run(self.target,self.source,self.log)
        self.assertTrue(reuse.pair_ready(self.target,'org1-org2'))
        self.assertTrue((self.target/sorted_name).exists())
        self.assertFalse((self.target/self.maf).exists())

    def test_other_run_log_cannot_certify_matching_filename(self):
        self.log.write_text(self.log.read_text().replace(str(self.source),str(self.root/'other')))
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse((self.target/self.maf).exists())

    def test_incomplete_training_block_is_rejected_even_with_matching_hash(self):
        (self.source/self.train).write_text('# substitution percent identity: 95.0\n')
        value=manifest.load(self.source)
        for item in value['artifacts']:
            if item['path']==self.train:item['sha256']=manifest.sha256(self.source/self.train)
        manifest.atomic_json(self.source/manifest.MANIFEST,value)
        reuse.import_run(self.target,self.source,self.log)
        self.assertFalse((self.target/self.maf).exists())

    def test_maf_coordinate_and_row_checks(self):
        for text in (MAF.replace('0 4 + 10','8 4 + 10'),MAF.replace('1 4 - 12','1 3 - 12'),
                     MAF.replace('ACGT\n\n','ACG\n\n'),MAF+'a\ns third 0 1 + 1 A\n'):
            with self.subTest(text=text):
                path=self.root/'bad.maf';path.write_text(text)
                with self.assertRaises(ValueError): reuse.validate_maf(path)


if __name__=='__main__': unittest.main()
