#!/usr/bin/env python3
"""Explicit, restartable PBS retry campaign. Original results are read-only inputs."""
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'src/lib'))
import run_manifest as runs

PILOTS = [('primates',6),('telluraves',16),('telluraves',5),('telluraves',18)]
TOOLS = ('lastal','lastdb','last-train','last-split','maf-sort','maf-join','maf-linked')


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.tmp.'+str(os.getpid()))
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)


def inventory(taxon, source):
    source = Path(source).resolve()
    settings = json.loads((source/'trio_selection/selection_run.json').read_text())
    with (source/'trio_selection/selected_trios.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    existing = {}
    for p in source.glob('*/*/metadata/metadata_manifest.json'):
        v = json.loads(p.read_text())
        key = tuple(o['accession'] for o in v['organisms'])
        if key in existing:
            raise ValueError('Ambiguous existing run: '+str(key))
        existing[key] = (str(p.parent.parent),v)
    result = []
    for n,row in enumerate(rows,1):
        accs = [row[k] for k in ('out_acc','in1_acc','in2_acc')]
        path,v = existing.get(tuple(accs),(None,{}))
        done = v.get('status') == 'complete'
        result.append(dict(key=f'{taxon}-{n}',taxon=taxon,trio=n,accessions=accs,
            source_job=str(source),source_run=path,settings=settings,
            date=settings['date'],threads=settings['threads'],
            status='existing_complete' if done else 'pending',
            pilot=(taxon,n) in PILOTS and not done,job_id=None,
            result_path=path if done else None))
    return result


def next_entries(entries, maximum):
    if any(e['status']=='submission_unknown' for e in entries):
        return []
    active = sum(e['status'] in ('submitted','running','submitting') for e in entries)
    capacity = max(0,maximum-active)
    pilots = [e for e in entries if e['pilot']]
    if any(e['status'] in ('failed','submission_unknown') for e in pilots):
        return []
    pool = [e for e in entries if e['status']=='pending']
    if any(e['status']!='complete' for e in pilots):
        pool = [e for e in pool if e['pilot']]
        pool.sort(key=lambda e:PILOTS.index((e['taxon'],e['trio'])))
    return pool[:capacity]


def parse_pbs(text):
    def field(name):
        m=re.search(r'^\s*'+re.escape(name)+r'\s*=\s*(.*?)\s*$',text,re.M)
        return m.group(1) if m else None
    state=field('job_state')
    if state is None:
        raise ValueError('PBS did not return a job state: '+text[-500:])
    rc=field('Exit_status')
    return dict(state=state,exit=int(rc) if rc is not None else None,
        memory=field('resources_used.mem'),walltime=field('resources_used.walltime'))


def query_pbs(job_id):
    # Miyabi -H selects history ONLY; it excludes queued/running jobs.
    responses=[]
    for flags in (['-f'],['-f','-H']):
        result=subprocess.run(['qstat',*flags,job_id],capture_output=True,text=True,timeout=60)
        responses.append(result.stdout+result.stderr)
        if result.returncode==0:
            try:return parse_pbs(result.stdout)
            except ValueError:pass
    raise RuntimeError('PBS state unavailable; no new submissions: '+' | '.join(responses))


def runtime_environment(campaign, entry):
    env=os.environ.copy()
    for key in ('EVOSUBSTER_REUSE_SOURCE','EVOSUBSTER_REUSE_LOG','EVOSUBSTER_METADATA_INPUT',
                'EVOSUBSTER_LAST_OPTIONS','EVOSUBSTER_SELECTION_SOURCE'):
        env.pop(key,None)
    env.update({str(k):str(v) for k,v in entry['settings']['last_environment'].items()})
    env.update(THREAD_NUM=str(entry['threads']),LASTAL_THREADS='1',LASTAL_PIPELINE_MODE='staged',
        OMP_NUM_THREADS=str(entry['threads']),OPENBLAS_NUM_THREADS='1',
        EVOSUBSTER_LAST_OPTIONS=json.dumps(entry['settings']['last_options']),
        EVOSUBSTER_SELECTION_SOURCE=str(Path(entry['source_job'])/'trio_selection/selected_trios.tsv'))
    env['PATH']=campaign['last_bin']+os.pathsep+env.get('PATH','')
    if entry['source_run']:
        env['EVOSUBSTER_REUSE_SOURCE']=entry['source_run']
        env['EVOSUBSTER_REUSE_LOG']=str(Path(entry['source_job'])/'job.log')
    tmp=Path(campaign['directory'])/'tmp'/entry['key']
    env.update(TMPDIR=str(tmp),TMP=str(tmp),TEMP=str(tmp))
    return env


def report(c):
    path=Path(c['directory'])/'summary.tsv'
    temp=path.with_suffix('.tmp')
    columns=['taxon','trio','status','job_id','last_version','source_run','result_path',
             'reused_pairs','gff_missing','memory','walltime','error']
    with temp.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,delimiter='\t',extrasaction='ignore')
        writer.writeheader()
        for e in c['entries']:
            r=dict(e)
            r['last_version']='1648' if e['status']=='existing_complete' else '1654'
            r.update({k:e.get('pbs',{}).get(k) for k in ('memory','walltime')})
            if e.get('result_path'):
                p=Path(e['result_path'])/'metadata/metadata_manifest.json'
                if p.exists():
                    v=json.loads(p.read_text())
                    r['reused_pairs']=','.join(sorted(v.get('legacy_reuse',{}).get('pairs',{})))
                    if r['reused_pairs']:r['last_version']='1654; inherited pairs=1648'
                    r['gff_missing']=','.join(o['slot'] for o in v['organisms'] if o.get('gff_path') is None)
            writer.writerow(r)
    temp.replace(path)


def save(c):
    c['updated_at']=datetime.now(timezone.utc).isoformat()
    atomic(Path(c['directory'])/'campaign.json',c)
    report(c)


def init(args):
    directory=Path(args.directory).resolve()
    if directory.exists():
        raise ValueError('Campaign directory already exists')
    repo=Path(args.repo).resolve()
    entries=inventory('primates',repo/'results/primates/3826407.opbs')
    entries+=inventory('telluraves',repo/'results/telluraves/3835636.opbs')
    expected={'primates':[1,3,6,8,9,10,11,12,13,14], 'telluraves':[5,16,18,19,20,21,22,23,24,25,26,27]}
    for taxon,nums in expected.items():
        found=[e['trio'] for e in entries if e['taxon']==taxon and e['status']=='pending']
        if found!=nums:raise ValueError(f'Retry inventory changed: {taxon}: {found}')
    validation={v['path']:v for v in json.loads(Path(args.validation).read_text())}
    for e in entries:
        if e['status']=='existing_complete':
            v=validation.get(e['source_run'],{})
            p=Path(e['source_run'])/'metadata/metadata_manifest.json'
            if v.get('status')!='validated' or runs.sha256(p)!=v.get('manifest_sha256'):
                raise ValueError('Completed source has not been validated: '+str(p))
    last_bin=repo/'.venv/tools/last-1654/bin'
    for name in TOOLS:
        if not (last_bin/name).is_file():raise ValueError('Missing LAST1654 tool: '+name)
    version=subprocess.check_output([str(last_bin/'lastal'),'--version'],text=True).strip()
    if version!='lastal 1654':raise ValueError('Unexpected LAST: '+version)
    tests=json.loads((last_bin.parent/'official-tests.json').read_text())
    required={'last-test.sh','last-train-test.sh','last-split-test.sh','maf-convert-test.sh','maf-cut-test.sh','maf-linked-test.sh','maf-swap-test.sh'}
    passed={t['test'] for t in tests if t['exit']==0}
    if not required.issubset(passed):raise ValueError('Required official LAST alignment/MAF tests must pass')
    directory.mkdir(parents=True)
    for name in ('logs','status','runs','tmp'): (directory/name).mkdir()
    code=directory/'code'
    for name in ('src','config','scripts'):
        shutil.copytree(repo/name,code/name,ignore=shutil.ignore_patterns('__pycache__'))
    tool_source=last_bin.parent
    last_bin=code/'tools/last-1654/bin'
    shutil.copytree(tool_source/'bin',last_bin)
    shutil.copyfile(tool_source/'source.json',last_bin.parent/'source.json')
    shutil.copyfile(tool_source/'official-tests.json',last_bin.parent/'official-tests.json')
    code_hashes={str(p.relative_to(code)):runs.sha256(p) for p in code.rglob('*') if p.is_file()}
    atomic(directory/'code-checksums.json',code_hashes)
    c=dict(schema=1,repo=str(repo),code_dir=str(code),directory=str(directory),last_bin=str(last_bin),
        max_active=4,phase='pilots',entries=entries,
        source=json.loads((last_bin.parent/'source.json').read_text()),
        created_at=datetime.now(timezone.utc).isoformat())
    shutil.copyfile(args.validation,directory/'original-validation.json')
    save(c)
    print(directory/'summary.tsv')


def run_entry(c,key):
    entry=next(e for e in c['entries'] if e['key']==key)
    if entry['status']=='existing_complete':raise ValueError('Completed originals cannot be rerun')
    root=Path(c['directory']); status_path=root/'status'/f'{key}.json'
    tmp=root/'tmp'/key;tmp.mkdir(exist_ok=True)
    lock=(tmp/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    env=runtime_environment(c,entry)
    result=root/'runs'/key
    log=root/'logs'/f'{key}.{os.environ.get("PBS_JOBID","local")}.log'
    state=dict(key=key,status='running',job_id=os.environ.get('PBS_JOBID'),started_at=datetime.now(timezone.utc).isoformat(),log=str(log))
    atomic(status_path,state)
    child=None
    def terminate(signum,frame):
        if child is not None:
            try:os.killpg(child.pid,signum)
            except ProcessLookupError:pass
        raise InterruptedError('Received signal '+str(signum))
    signal.signal(signal.SIGTERM,terminate)
    signal.signal(signal.SIGINT,terminate)
    try:
        recorded=json.loads((root/'code-checksums.json').read_text())
        for relative,digest in recorded.items():
            if runs.sha256(Path(c['code_dir'])/relative)!=digest:
                raise ValueError('Execution code changed: '+relative)
        for tool in TOOLS:
            found=shutil.which(tool,path=env['PATH'])
            if found!=str(Path(c['last_bin'])/tool):raise ValueError('Mixed LAST tool: '+tool)
        command=['bash',str(Path(c['code_dir'])/'src/sbst_fromDwl.sh'),entry['date'],*entry['accessions'],
                 '--genome-dir',entry['settings']['genome_dir'],'--out-dir',str(result),
                 '--thread',str(entry['threads']),'--train-cache-dir',str(root/'cache'/key)]
        with log.open('a',buffering=1) as stream:
            stream.write(json.dumps(dict(command=command,last_version='1654',tools={t:shutil.which(t,path=env['PATH']) for t in TOOLS},
                                    lastal_threads=1,training_threads=entry['threads'],mode='staged'))+'\n')
            stream.flush()
            child=subprocess.Popen(command,cwd=c['repo'],env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            rc=child.wait()
        if rc:raise RuntimeError('pipeline exit='+str(rc))
        paths=list(result.glob('*/*/metadata/metadata_manifest.json'))
        if len(paths)!=1:raise ValueError('Expected exactly one result manifest')
        run=paths[0].parent.parent
        runs.validated(run,complete=True)
        state.update(status='complete',exit=0,result_path=str(run),validated=True)
    except BaseException as error:
        state.update(status='failed',exit=1,error=str(error))
    state['finished_at']=datetime.now(timezone.utc).isoformat()
    atomic(status_path,state)
    return state['exit']


def advance(c):
    root=Path(c['directory'])
    for e in c['entries']:
        if e['status']=='submitting':
            e.update(status='submission_unknown',error='Submission interrupted; reconcile PBS before retrying')
        if e['status'] not in ('submitted','running'):continue
        e['pbs']=query_pbs(e['job_id'])
        if e['pbs']['state']!='F':
            e['status']='running' if e['pbs']['state']=='R' else 'submitted'
            continue
        state_path=root/'status'/f"{e['key']}.json"
        state=json.loads(state_path.read_text()) if state_path.exists() else {}
        if e['pbs']['exit']==0 and state.get('status')=='complete' and state.get('validated'):
            e.update(status='complete',result_path=state['result_path'])
        else:
            e.update(status='failed',error=state.get('error',f"PBS exit={e['pbs']['exit']}"))
        save(c)
    pilots=[e for e in c['entries'] if e['pilot']]
    if any(e['status'] in ('failed','submission_unknown') for e in pilots):c['phase']='pilot_failed'
    elif all(e['status']=='complete' for e in pilots):c['phase']='remaining'
    for e in next_entries(c['entries'],c['max_active']):
        e['status']='submitting';save(c)
        cmd=['qsub','-N','retry_'+e['key'].replace('-','_'),'-v',
             'CAMPAIGN='+c['directory']+',TRIO_KEY='+e['key']+',PROJECT_DIR='+c['repo'],
             '-o',str(root/'logs'/f"{e['key']}.pbs.log"),
             str(Path(c['code_dir'])/'scripts/retry_trio_miyabi_g.pbs')]
        try:
            env=os.environ.copy();env['TMPDIR']='/tmp'
            p=subprocess.run(cmd,cwd=c['repo'],env=env,capture_output=True,text=True,timeout=60)
            if p.returncode:
                e.update(status='submission_unknown',error='qsub response requires reconciliation: '+p.stdout+p.stderr)
            else:
                match=re.search(r'(?m)^\s*(\d+(?:\.[A-Za-z0-9_.-]+)?)\s*$',p.stdout)
                if not match:raise ValueError('Unrecognized qsub response: '+p.stdout)
                e.update(status='submitted',job_id=match.group(1))
                print('submitted',e['key'],e['job_id'],flush=True)
        except Exception as error:
            e.update(status='submission_unknown',error=str(error));save(c);raise
        save(c)
        if e['status']!='submitted':break
    if all(e['status'] in ('existing_complete','complete','failed') for e in c['entries']):
        c['phase']='finished' if all(e['status']!='failed' for e in c['entries']) else 'finished_with_failures'
    save(c)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('init');p.add_argument('directory');p.add_argument('--repo',default=str(REPO));p.add_argument('--validation',required=True)
    for name in ('advance','watch','report','run'):
        p=sub.add_parser(name);p.add_argument('directory')
        if name=='run':p.add_argument('key')
    args=parser.parse_args()
    if args.command=='init':init(args);return
    directory=Path(args.directory).resolve()
    if args.command=='run':return run_entry(json.loads((directory/'campaign.json').read_text()),args.key)
    lock=(directory/'controller.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    while True:
        c=json.loads((directory/'campaign.json').read_text())
        if args.command=='report':report(c);return
        try:advance(c)
        except Exception as error:
            print(datetime.now(timezone.utc).isoformat(),str(error),file=sys.stderr,flush=True)
            if args.command!='watch':raise
        if args.command!='watch':return
        if c['phase'].startswith('finished'):return
        if c['phase']=='pilot_failed' and not any(e['status'] in ('submitted','running') for e in c['entries']):return 1
        time.sleep(60)


if __name__=='__main__':
    sys.exit(main() or 0)
