#!/usr/bin/env python3
"""Read-only CLI preview. Never downloads, allocates runs, trains, or publishes files."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import run_manifest as manifest
import selection_run

ROOT=Path(__file__).resolve().parents[2]


def root_path(value):
    path=Path(value)
    return path.resolve() if path.is_absolute() else (ROOT/path).resolve()


def short_name(name,slot):
    parts=name.split('_')
    first=parts[0]
    second=parts[1] if len(parts)>1 else first[3:]
    return first[:3]+second[:3]+str(slot)


def defaults():
    text=subprocess.check_output(['yq','eval','-o=json','.',str(ROOT/'config/dwl_config.yaml')],text=True)
    return json.loads(text)['paths']


def preview_local(date,fastas,gffs,accessions,out_dir,threads,extra=None):
    missing=[str(p) for p in fastas if not Path(p).is_file()]
    missing += [str(p) for p in gffs if p!='NO_GFF_FILE' and not Path(p).is_file()]
    if missing:
        print(json.dumps({'action':'unavailable','reason':'Input not acquired or unreadable; no downloads performed','inputs':missing}))
        return
    os.environ['THREAD_NUM']=str(threads)
    inputs=[]
    for slot,(fasta,gff,acc) in enumerate(zip(fastas,gffs,accessions),1):
        inputs += [str(Path(fasta).resolve()), str(Path(gff).resolve()) if gff!='NO_GFF_FILE' else gff,
                   acc or '', short_name(Path(fasta).parent.name,slot)]
    args=argparse.Namespace(label=date,threads=int(threads),inputs=inputs)
    value=manifest.descriptor(args)
    if extra:
        for entry,metadata in zip(value['organisms'],extra):
            entry.update(metadata)
    base=root_path(out_dir or defaults()['out_dir'])
    existing=manifest.resume_matches(base,value)
    source=os.environ.get('EVOSUBSTER_REUSE_FROM') or None
    if existing:
        previous=manifest.load(existing)
        value['stage_records']=previous.get('stage_records',{})
        for current,old in zip(value['organisms'],previous['organisms']):
            current['short_name']=old['short_name']
    import stage_reuse
    result=stage_reuse.plan(value,source,os.environ.get('EVOSUBSTER_REUSE_LAST_VERSION','same'),target_run=existing)
    if existing:
        old=manifest.load(existing)
        if old['status']=='complete':
            try:
                manifest.validated(existing,complete=True)
                valid=all(stage_reuse.ready(existing,r['kind'],r['slot']) for r in old.get('stage_records',{}).values())
                result['run_action']='skip_complete' if valid else 'resume_incomplete'
            except (ValueError,OSError,KeyError):
                result['run_action']='resume_incomplete'
    result.update(date=date,accessions=accessions,output_root=str(base),existing_run=str(existing) if existing else None)
    for row in result['stages']:
        record=row.pop('record',None)
        if record:
            row['producer_last_version']=record['producer'].get('last_version')
    print(json.dumps(result,indent=2))


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out-dir',default='')
    p.add_argument('--genome-dir','--base-genomes',dest='genome_dir',default='')
    p.add_argument('--train-cache-dir',default='')
    p.add_argument('--thread',default='')
    p.add_argument('--org2-gff',default='NO_GFF_FILE')
    p.add_argument('--org3-gff',default='NO_GFF_FILE')
    p.add_argument('--accession-ids',nargs=3)
    p.add_argument('--trios',default='')
    p.add_argument('--tree',default='')
    p.add_argument('--selection-config',default='')
    p.add_argument('--idt-only',action='store_true')
    p.add_argument('--force',action='store_true')
    p.add_argument('positional',nargs='*')
    return p


def download_preview(args):
    date=args.positional[0] if args.positional else ''
    trios=[]
    if args.tree and args.trios:
        raise ValueError('--tree and --trios cannot be combined')
    if (args.tree or args.trios) and len(args.positional)>1:
        raise ValueError('Tree/trios mode does not accept positional accessions')
    if args.tree:
        print(json.dumps({'action':'unavailable','reason':'Tree selection must be resolved first; preview never runs selection. Use --trios with a saved selected_trios.tsv to preview exact trios.'}))
        return
    if args.trios:
        trios=selection_run.read_trios(args.trios)
        resolve=[sys.executable,str(ROOT/'src/lib/selection_run.py'),'resolve','--trios',args.trios,
                 '--config',args.selection_config,'--date',date,'--threads',args.thread,
                 '--genome-dir',args.genome_dir,'--out-dir',args.out_dir,'--train-cache-dir',args.train_cache_dir]
        for tool in os.environ.get('EVOSUBSTER_EXPLICIT_LAST','').split():
            resolve += ['--explicit-last',tool]
        fields=subprocess.check_output(resolve).decode().rstrip('\0').split('\0')
        date,args.genome_dir,args.out_dir,_,args.thread=fields[:5]
        os.environ.update(EVOSUBSTER_LAST_OPTIONS=fields[5],LASTDB_ALIGNMENT_SEED=fields[6],
                          LASTAL_QUERY_BATCH_SIZE=fields[7],LASTAL_SPLIT_MEMORY=fields[8])
    else:
        if len(args.positional)!=4:
            raise ValueError('Expected DATE and three accessions')
        trios=[args.positional[1:]]
    genomes=root_path(args.genome_dir or defaults()['base_genomes'])
    bundles={}
    for path in genomes.rglob('assembly_bundle.json'):
        value=json.loads(path.read_text())
        accession=value.get('accession')
        if accession:
            bundles.setdefault(accession,[]).append(path)
    for trio in trios:
        fastas=[];gffs=[]
        unavailable=[]
        for accession in trio:
            candidates=bundles.get(accession,[])
            if len(candidates)>1:
                raise ValueError('Ambiguous local assembly bundles: '+accession)
            if not candidates:
                unavailable.append(accession)
                continue
            path=candidates[0]
            data=json.loads(path.read_text())
            fasta=manifest.artifact_path(path.parent,data['fasta'])
            fastas.append(str(fasta))
            gffs.append(manifest.bundle_gff(fasta,accession))
        if unavailable:
            print(json.dumps({'accessions':trio,'action':'unavailable',
                              'reason':'Inputs not acquired; no downloads performed','inputs':unavailable}))
            continue
        preview_local(date,fastas,gffs,trio,args.out_dir,args.thread or '8')


def main():
    entry=sys.argv[1]
    args,unknown=parser().parse_known_args(sys.argv[2:])
    if unknown:
        raise ValueError('Unsupported preview arguments: '+' '.join(unknown))
    if entry=='from-download':
        download_preview(args)
    else:
        if len(args.positional)!=5:
            raise ValueError('Expected DATE, three FASTAs and outgroup GFF')
        date,*values=args.positional
        fastas=values[:3]
        accessions=args.accession_ids
        if not accessions:
            import re
            accessions=[]
            for fasta in fastas:
                found=re.match(r'(GC[AF]_[0-9]+\.[0-9]+)(?:_|\.|$)',Path(fasta).name)
                accessions.append(found.group(1) if found else '')
        preview_local(date,fastas,[values[3],args.org2_gff,args.org3_gff],accessions,args.out_dir,args.thread or '8')

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError,subprocess.CalledProcessError) as error:
        print('Error: '+str(error),file=sys.stderr)
        sys.exit(1)
