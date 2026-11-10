"""Content-based stage reuse must never certify an incomplete or changed stage."""
import copy
import fcntl
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
import run_manifest as manifest
try:
    import stage_reuse as reuse
except ImportError:
    reuse = None

TRAIN = '''# substitution percent identity: 95.0
# substitution percent identity: 99.0
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


class StageReuseTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(reuse, 'stage_reuse module is required')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source, self.target = self.root/'source', self.root/'target'
        organisms = []
        for i in range(1, 4):
            fasta = self.root/f'genome{i}.fa'; fasta.write_text(f'>seq{i}\nACGT\n')
            gff = self.root/f'genome{i}.gff'; gff.write_text('seq\tx\tCDS\t1\t4\t.\t+\t0\tID=x\n')
            organisms.append(dict(slot=f'org{i}', accession=f'any{i}', short_name=f'Sp{i}',
                fasta_path=str(fasta), fasta_sha256=manifest.sha256(fasta),
                gff_path=str(gff), gff_sha256=manifest.sha256(gff)))
        last = dict(schema=1,lastdb_train=['-P','8','-u','RY4'],last_train=['-P','8','--revsym'],
            lastdb_align=['-P','8','-u','RY128'],lastal=['-P','8','-C','2'],last_split=['-r'])
        self.value = dict(schema_version=2, run_key='source-key',trio_key='trio',date='old',created_at='now',
            organisms=organisms,artifacts=[],status='pending',options=dict(last=last,last_version='lastal 1654',
            algorithm='pipeline-v2',config={'names':'ignored'}))
        manifest.prepare(self.source,copy.deepcopy(self.value))
        target = copy.deepcopy(self.value); target['run_key']='target-key'; target['date']='new'
        for o in target['organisms']:o['short_name']='Renamed'+o['slot']
        manifest.prepare(self.target,target)

    def write(self, kind, slot, outputs):
        result = {}
        for role, data in outputs.items():
            path = f'intermediateFiles/{kind}-{slot}.{role}'
            dest = self.source/path; dest.parent.mkdir(exist_ok=True); dest.write_text(data)
            result[role]=path
        reuse.register(self.source,kind,slot,result)
        return result

    def stage(self, value, kind, slot):
        return next(s for s in value['stages'] if s['kind']==kind and s['slot']==slot)

    def pair(self):
        train=self.write('train','org1-org2',{'train':TRAIN})
        return self.write('pair','org1-org2',{'maf':MAF,'train':TRAIN}), train

    def test_explicit_completion_and_hash_verification(self):
        self.assertFalse(reuse.ready(self.source,'gc','org1'))
        names=self.write('gc','org1',{'gc':'0.5\n'})
        self.assertTrue(reuse.ready(self.source,'gc','org1'))
        (self.source/names['gc']).write_text('changed')
        self.assertFalse(reuse.ready(self.source,'gc','org1'))

    def test_registered_and_reused_training_keeps_penultimate_provenance(self):
        self.write('train', 'org1-org2', {'train': TRAIN})
        source_value = manifest.validated(self.source)
        artifact = next(a for a in source_value['artifacts'] if a['kind'] == 'train')
        self.assertEqual(artifact['identity_percent'], 95.0)
        self.assertEqual(artifact['identity_source'], 'last-train-penultimate')
        reuse.apply(self.target, self.source)
        target_value = manifest.validated(self.target)
        reused = next(a for a in target_value['artifacts'] if a['kind'] == 'train')
        self.assertEqual(reused['identity_percent'], 95.0)
        self.assertEqual(reused['identity_source'], 'last-train-penultimate')

    def test_apply_maps_names_and_preserves_source_and_original_producer(self):
        self.pair(); self.write('gc','org1',{'gc':'0.5\n'})
        before={str(p.relative_to(self.source)):p.read_bytes() for p in self.source.rglob('*') if p.is_file()}
        reuse.apply(self.target,self.source)
        self.assertTrue(reuse.ready(self.target,'pair','org1-org2'))
        self.assertTrue(reuse.ready(self.target,'train','org1-org2'))
        maf=self.target/'intermediateFiles/Renamedorg12Renamedorg2_one2one_new.maf'
        self.assertEqual(maf.read_text(),MAF)
        self.assertTrue((self.target/'statistics/misc/Renamedorg1_gcContent_new.out').exists())
        self.assertEqual(before,{str(p.relative_to(self.source)):p.read_bytes() for p in self.source.rglob('*') if p.is_file()})
        self.assertNotEqual(maf.stat().st_ino,(self.source/'intermediateFiles/pair-org1-org2.maf').stat().st_ino)
        record=manifest.load(self.target)['stage_records']['pair:org1-org2']
        self.assertEqual(record['producer']['last_version'],'lastal 1654')
        self.assertEqual(record['provenance']['source_run'],str(self.source))
        manifest.validated(self.target)

    def test_scoped_input_gff_and_seed_invalidation(self):
        self.pair(); self.write('gc','org1',{'gc':'0.5'}); self.write('gc','org3',{'gc':'0.5'})
        self.write('coding','org1',{'coding':'100'})
        value=manifest.load(self.target)
        value['organisms'][0]['gff_sha256']='changed'
        value['organisms'][2]['fasta_sha256']='changed'
        value['options']['last']['lastdb_align'][-1]='YASS'
        plan=reuse.plan(value,self.source)
        for kind,slot,action in [('gc','org1','reuse'),('gc','org3','recompute'),
            ('coding','org1','recompute'),('train','org1-org2','reuse'),('pair','org1-org2','recompute')]:
            self.assertEqual(self.stage(plan,kind,slot)['action'],action)

    def test_version_policy_imports_pair_bundle_but_not_independent_training(self):
        self.pair(); self.write('train','org2-org3',{'train':TRAIN})
        value=manifest.load(self.target); value['options']['last_version']='lastal 9999'
        manifest.atomic_json(self.target/manifest.MANIFEST,value)
        self.assertEqual(self.stage(reuse.plan(value,self.source),'pair','org1-org2')['action'],'recompute')
        plan=reuse.apply(self.target,self.source,'allow-existing-pairs')
        self.assertEqual(self.stage(plan,'pair','org1-org2')['action'],'reuse')
        self.assertEqual(self.stage(plan,'train','org2-org3')['action'],'recompute')
        self.assertTrue(reuse.ready(self.target,'pair','org1-org2'))
        self.assertTrue(reuse.ready(self.target,'train','org1-org2'))
        self.assertEqual(manifest.load(self.target)['stage_records']['pair:org1-org2']['producer']['last_version'],'lastal 1654')

    def test_lastal_threads_ignored_but_training_threads_strict(self):
        self.pair()
        value=manifest.load(self.target); value['options']['last']['lastal'][1]='1'
        self.assertEqual(self.stage(reuse.plan(value,self.source),'pair','org1-org2')['action'],'reuse')
        value['options']['last']['last_train'][1]='1'
        self.assertEqual(self.stage(reuse.plan(value,self.source),'pair','org1-org2')['action'],'recompute')

    def test_ambiguity_rejected_and_corruption_recomputed(self):
        outputs=self.write('gc','org1',{'gc':'0.5'})
        shutil.copytree(self.source,self.root/'second')
        with self.assertRaisesRegex(ValueError,'[Aa]mbigu'):reuse.plan(manifest.load(self.target),self.root,target_run=self.target)
        shutil.rmtree(self.root/'second')
        (self.source/outputs['gc']).write_text('bad')
        item=self.stage(reuse.plan(manifest.load(self.target),self.source),'gc','org1')
        self.assertEqual(item['action'],'recompute')
        self.assertTrue(item['reason'])

    def test_active_explicit_source_rejected_without_creating_lock(self):
        self.write('gc','org1',{'gc':'0.5'})
        reuse.plan(manifest.load(self.target),self.source)
        self.assertFalse((self.source/'.run.lock').exists())
        with (self.source/'.run.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with self.assertRaisesRegex(ValueError,'active|locked'):reuse.plan(manifest.load(self.target),self.source)

    def test_reject_partial_maf_bad_training_and_escaped_outputs(self):
        self.write('train','org1-org2',{'train':TRAIN})
        path=self.source/'bad.maf';path.write_text(MAF+'a\ns x 0 1 + 2 A\n')
        train=self.source/'ok.train';train.write_text(TRAIN)
        with self.assertRaises(ValueError):reuse.register(self.source,'pair','org1-org2',{'maf':'bad.maf','train':'ok.train'})
        train.write_text('# substitution percent identity: 95.0\n')
        with self.assertRaises(ValueError):reuse.register(self.source,'train','org1-org2',{'train':'ok.train'})
        with self.assertRaises(ValueError):reuse.register(self.source,'gc','org1',{'gc':'../genome1.fa'})

    def test_no_overwrite_or_symlink_source_discovery(self):
        self.write('gc','org1',{'gc':'0.5'})
        destination=self.target/'statistics/misc/Renamedorg1_gcContent_new.out'
        destination.parent.mkdir(parents=True);destination.write_text('other')
        with self.assertRaisesRegex(ValueError,'overwrite'):reuse.apply(self.target,self.source)
        other=self.root/'scan';other.mkdir();(other/'external').symlink_to(self.source,target_is_directory=True)
        self.assertEqual(self.stage(reuse.plan(manifest.load(self.target),other),'gc','org1')['action'],'recompute')

    def test_relocation_keeps_pair_ready_and_imports_sorted_target_name(self):
        outputs, _ = self.pair()
        old=outputs['maf']; new=old[:-4]+'_sorted.maf'
        (self.source/old).rename(self.source/new)
        self.assertTrue(reuse.relocate(self.source,old,new))
        self.assertTrue(reuse.ready(self.source,'pair','org1-org2'))
        reuse.apply(self.target,self.source)
        self.assertTrue((self.target/'intermediateFiles/Renamedorg12Renamedorg2_one2one_new_sorted.maf').exists())
        self.assertTrue(reuse.ready(self.target,'pair','org1-org2'))

    def test_repeat_register_preserves_imported_old_producer(self):
        self.pair()
        value=manifest.load(self.target); value['options']['last_version']='lastal 9999'
        manifest.atomic_json(self.target/manifest.MANIFEST,value)
        reuse.apply(self.target,self.source,'allow-existing-pairs')
        record=manifest.load(self.target)['stage_records']['pair:org1-org2']
        reuse.register(self.target,'pair','org1-org2',{o['role']:o['path'] for o in record['outputs']})
        self.assertEqual(record,manifest.load(self.target)['stage_records']['pair:org1-org2'])

    def test_source_change_during_copy_aborts_before_outputs_are_published(self):
        self.write('gc','org1',{'gc':'0.5'})
        original=shutil.copyfile
        def changed(source,dest):
            result=original(source,dest)
            value=manifest.load(self.source); value['status']='identity_complete'
            manifest.atomic_json(self.source/manifest.MANIFEST,value)
            return result
        with patch.object(reuse.shutil,'copyfile',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'[Ss]ource manifest changed'):reuse.apply(self.target,self.source)
        self.assertFalse((self.target/'statistics/misc/Renamedorg1_gcContent_new.out').exists())
        self.assertFalse(manifest.load(self.target).get('stage_records'))

    def test_manifest_symlink_cannot_escape_source_root(self):
        self.write('gc','org1',{'gc':'0.5'})
        malicious=self.root/'malicious'; (malicious/'metadata').mkdir(parents=True)
        (malicious/manifest.MANIFEST).symlink_to(self.source/manifest.MANIFEST)
        with self.assertRaisesRegex(ValueError,'[Ss]ymlink'):reuse.plan(manifest.load(self.target),malicious)

    def test_verified_native_seals_survive_manifest_preparation(self):
        path=self.source/'intermediateFiles/profile.train'; path.parent.mkdir(); path.write_text(TRAIN)
        manifest.seal_file(path)
        reuse.register(self.source,'train','org1-org2',{'train':str(path.relative_to(self.source))})
        manifest.prepare(self.source,manifest.load(self.source))
        self.assertTrue(Path(str(path)+'.complete.json').exists())
        self.assertTrue(reuse.ready(self.source,'train','org1-org2'))

    def test_same_run_identical_genomes_do_not_make_metric_reuse_ambiguous(self):
        value=manifest.load(self.source)
        value['organisms'][1]['fasta_sha256']=value['organisms'][0]['fasta_sha256']
        manifest.atomic_json(self.source/manifest.MANIFEST,value)
        self.write('gc','org1',{'gc':'0.5'})
        self.write('gc','org2',{'gc':'0.5'})
        plan=reuse.plan(manifest.load(self.target),self.source)
        self.assertEqual(self.stage(plan,'gc','org1')['action'],'reuse')
        self.assertEqual(self.stage(plan,'gc','org1')['source_slot'],'org1')

    def test_source_metadata_corruption_is_rejected(self):
        self.write('gc','org1',{'gc':'0.5'})
        (self.source/'metadata/org1_assembly.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'Metadata checksum'):reuse.plan(manifest.load(self.target),self.source)

    def test_imported_training_gets_a_valid_seal(self):
        self.write('train','org1-org2',{'train':TRAIN})
        reuse.apply(self.target,self.source)
        path=self.target/'intermediateFiles/Renamedorg12Renamedorg2_new.train'
        manifest.seal_file(path,verify=True)
        self.assertIn(str(path.relative_to(self.target))+'.complete.json',[a['path'] for a in manifest.load(self.target)['artifacts']])

    def test_normalized_alignment_guard_is_preserved(self):
        self.write('train','org1-org2',{'train':TRAIN})
        maf=self.source/'intermediateFiles/native.maf'; maf.write_text(MAF)
        train=self.source/'intermediateFiles/train-org1-org2.train'
        from last_options import alignment_settings
        guard=Path(str(maf)+'.last_options.json')
        guard.write_text(json.dumps(alignment_settings(manifest.load(self.source)['options']['last'])))
        reuse.register(self.source,'pair','org1-org2',{'maf':str(maf.relative_to(self.source)),'train':str(train.relative_to(self.source))})
        manifest.prepare(self.source,manifest.load(self.source))
        self.assertTrue(guard.exists())

    def test_target_metadata_symlink_is_rejected_before_import(self):
        self.write('gc','org1',{'gc':'0.5'})
        metadata=self.target/'metadata'; moved=self.root/'outside-metadata'
        metadata.rename(moved); metadata.symlink_to(moved,target_is_directory=True)
        before=(moved/'metadata_manifest.json').read_bytes()
        with self.assertRaisesRegex(ValueError,'[Ss]ymlink'):reuse.apply(self.target,self.source)
        self.assertEqual(before,(moved/'metadata_manifest.json').read_bytes())

    def test_changed_seed_definition_file_invalidates_pair(self):
        seed=self.root/'seed.txt'
        seed.write_text('1110101\n')
        value=manifest.load(self.source)
        value['options']['last']['lastdb_align']=['-P','8','-u',str(seed)]
        manifest.atomic_json(self.source/manifest.MANIFEST,value)
        self.pair()
        self.assertTrue(reuse.ready(self.source,'pair','org1-org2'))
        seed.write_text('1111111\n')
        self.assertFalse(reuse.ready(self.source,'pair','org1-org2'))


if __name__=='__main__':unittest.main()
