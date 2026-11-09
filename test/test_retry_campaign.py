import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('retry_campaign', ROOT/'scripts/retry_campaign.py')


class RetryCampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def module(self):
        self.assertTrue(Path(SPEC.origin).exists(), 'retry campaign controller is required')
        mod = importlib.util.module_from_spec(SPEC)
        SPEC.loader.exec_module(mod)
        return mod

    def test_inventory_excludes_completed_runs_and_marks_pilots(self):
        mod = self.module()
        source = self.root/'source'
        sel = source/'trio_selection'
        sel.mkdir(parents=True)
        rows = ['out_acc\tin1_acc\tin2_acc']
        for n in range(1,8):
            rows.append(f'GCF_{n}1.1\tGCF_{n}2.1\tGCF_{n}3.1')
        (sel/'selected_trios.tsv').write_text('\n'.join(rows)+'\n')
        (sel/'selection_run.json').write_text(json.dumps({'date':'20261108','threads':48,'genome_dir':str(self.root/'genomes'),'last_options':{},'last_environment':{}}))
        p=source/'key'/'date-key'/'metadata'
        p.mkdir(parents=True)
        (p/'metadata_manifest.json').write_text(json.dumps({'status':'complete','organisms':[{'accession':f'GCF_2{i}.1'} for i in range(1,4)]}))
        entries=mod.inventory('primates',source)
        self.assertEqual([e['trio'] for e in entries if e['status']=='existing_complete'],[2])
        self.assertEqual([e['trio'] for e in entries if e['pilot']],[6])
        self.assertEqual(len(entries),7)

    def test_pilots_gate_remaining_and_capacity_is_bounded(self):
        mod=self.module()
        es=[{'key':str(i),'status':'pending','pilot':i<4, 'taxon':mod.PILOTS[i][0] if i<4 else 'primates', 'trio':mod.PILOTS[i][1] if i<4 else i+20} for i in range(6)]
        self.assertEqual(mod.next_entries(es,4),es[:4])
        es[0]['status']='submitted'
        self.assertEqual(len(mod.next_entries(es,4)),3)
        for e in es[:4]:e['status']='complete'
        self.assertEqual(mod.next_entries(es,4),es[4:])
        es[0]['status']='failed'
        self.assertEqual(mod.next_entries(es,4),[])

    def test_unknown_submission_blocks_new_jobs(self):
        mod=self.module()
        entries=[{'status':'submission_unknown','pilot':False},{'status':'pending','pilot':False}]
        self.assertEqual(mod.next_entries(entries,4),[])

    def test_nonzero_qsub_is_uncertain_not_safe_to_resubmit(self):
        mod=self.module()
        c={'directory':str(self.root),'repo':str(self.root),'code_dir':str(self.root),'max_active':4,'phase':'pilots',
           'entries':[{'key':'primates-6','taxon':'primates','trio':6,'pilot':True,'status':'pending'}]}
        with patch.object(mod.subprocess,'run',return_value=subprocess.CompletedProcess([],1,'accepted 123.opbs','connection lost')):
            mod.advance(c)
        self.assertEqual(c['entries'][0]['status'],'submission_unknown')
        self.assertEqual(mod.next_entries(c['entries'],4),[])

    def test_live_query_precedes_history_and_falls_back_after_finish(self):
        mod=self.module()
        self.assertTrue(hasattr(mod,'query_pbs'), 'Live/history fallback required')
        with patch.object(mod.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'job_state = Q','')) as call:
            self.assertEqual(mod.query_pbs('123.opbs')['state'],'Q')
            self.assertNotIn('-H',call.call_args[0][0])
            self.assertEqual(call.call_args.kwargs.get('env',{}).get('PATH'),'/usr/local/bin:/opt/pbs/bin:/usr/bin:/bin')
        replies=[subprocess.CompletedProcess([],0,'No matching job found.',''),subprocess.CompletedProcess([],0,'job_state = F\nExit_status = 0','')]
        with patch.object(mod.subprocess,'run',side_effect=replies) as call:
            self.assertEqual(mod.query_pbs('123.opbs')['exit'],0)
            self.assertIn('-H',call.call_args[0][0])

    def test_pbs_requires_explicit_terminal_success(self):
        mod=self.module()
        self.assertEqual(mod.parse_pbs('job_state = F\nExit_status = -29\nresources_used.mem = 100kb')['exit'],-29)
        self.assertEqual(mod.parse_pbs('job_state = R\n')['state'],'R')
        with self.assertRaises(ValueError):mod.parse_pbs('Unknown Job Id')

    def test_runtime_uses_new_last_and_preserves_training_threads(self):
        mod=self.module()
        entry={'threads':48,'settings':{'last_environment':{},'last_options':{}},'source_run':None,'source_job':'/old','key':'primates-6'}
        env=mod.runtime_environment({'repo':str(self.root),'last_bin':str(self.root/'last/bin'),'directory':str(self.root/'campaign')},entry)
        self.assertEqual(env['LASTAL_THREADS'],'1')
        self.assertEqual(env['LASTAL_PIPELINE_MODE'],'staged')
        self.assertEqual(env['THREAD_NUM'],'48')
        self.assertTrue(env['PATH'].startswith(str(self.root/'last/bin')+':'))


if __name__=='__main__':unittest.main()
