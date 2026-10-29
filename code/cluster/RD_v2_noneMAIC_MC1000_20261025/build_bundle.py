"""Package eight RD v2 none+MAIC cells using the unchanged final scientific runtime."""
from pathlib import Path
import argparse, csv, hashlib, json, shutil, tarfile
P=Path(__file__).resolve().parent
BASE=P.parent/'RD_v2_full24_MC1000_20261024'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def replace_once(s,a,b):
    assert s.count(a)==1,(a[:90],s.count(a));return s.replace(a,b,1)

def stage():
    assert not (P/'VERSION.txt').exists(),'Already staged; use --package.'
    runtime=[]
    for p in sorted((BASE/'code').rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE/'code');q=P/'code'/rel
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        runtime.append(dict(path=rel.as_posix(),md5=hashlib.md5(q.read_bytes()).hexdigest(),sha256=sha(q)))
    table(P/'RUNTIME_FILES.csv',runtime)
    shutil.copyfile(BASE/'scenarios.csv',P/'scenarios.csv')
    cases=list(csv.DictReader((P/'scenarios.csv').open(encoding='utf-8')))
    rows=[]
    for c in cases:
        rows.append(dict(task_id=c['scenario_id'],**c,analysis='none_MAIC',estimation_mode='none',
            model_weighting='maic',ess_source='diffuse',em_level='high',weight_mode='all',nsim=1000,first_seed=1,
            cell=f"S{int(c['scenario_id']):02d}_{c['event']}_n{c['n']}_ess{round(100*float(c['ess_ratio'])):03d}_none_MAIC"))
    table(P/'tasks.csv',rows)
    table(P/'resources.csv',[dict(analysis='none_MAIC',n=n,cpus=64,memory='10G',walltime='03:00:00')for n in(100,500)])
    for name in ('common.R','prepare.R','run.R','collect.R','submit.sh','worker.sh','collect.sh'):
        s=(BASE/name).read_text('utf-8').replace('RDFULL','RDNM')
        s=s.replace('1:24','1:8').replace('1..24','1..8').replace('ALL 24 CELLS','ALL 8 CELLS')
        s=s.replace('RDv2_1000_g','RDnoneMAIC_g').replace('RDv2_report','RDnoneMAIC_report').replace('RDfull_','RDnoneMAIC_')
        if name=='common.R':
            start=s.index('tasks<-function(home){');end=s.index('settings.for<-function(',start)
            s=s[:start]+'''tasks<-function(home){
    x<-read.csv(file.path(home,'tasks.csv'),stringsAsFactors=FALSE)
    scenarios<-read.csv(file.path(home,'scenarios.csv'),stringsAsFactors=FALSE)
    must(nrow(x)==8L&&identical(x$task_id,1:8)&&identical(x$scenario_id,1:8)&&
        identical(x[c('scenario_id','event','n','ess_ratio')],scenarios)&&
        all(x$analysis=='none_MAIC')&&all(x$estimation_mode=='none')&&all(x$model_weighting=='maic')&&
        all(x$nsim==1000L)&&all(x$first_seed==1L)&&all(x$ess_source=='diffuse')&&
        all(x$em_level=='high')&&all(x$weight_mode=='all')&&!anyDuplicated(x$cell),
        'Invalid eight-case none+MAIC design')
    x
}
'''+s[end:]
            start=s.index('    prefix<-switch(r$analysis,');end=s.index('    must(is.finite(z$cpus)',start)
            s=s[:start]+'''    z$cpus<-as.integer(Sys.getenv('RDNM_CPUS',as.character(z$cpus)))
    z$memory<-Sys.getenv('RDNM_MEM',z$memory)
    z$walltime<-Sys.getenv('RDNM_TIME',z$walltime)
'''+s[end:]
        if name=='run.R':
            s=replace_once(s,"'1000 seeds;',cores,'workers; RD v2; final CI 600/1000/1000\\n'",
                "'1000 seeds;',cores,'workers; RD v2 none+MAIC; final CI 600/1000/1000\\n'")
        (P/name).write_text(s,encoding='utf-8',newline='\n')
    (P/'env.sh').write_text('''#!/bin/bash
export RDNM_HOME
RDNM_HOME=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export RDNM_OUT="${RDNM_OUT:-${SCRATCH:-/scratch/$USER}/RD_v2_noneMAIC_MC1000_results_20261025}"
export RDNM_ACCOUNT="${RDNM_ACCOUNT:-def-liteep}"
# Optional overrides: RDNM_CPUS, RDNM_MEM, RDNM_TIME. Defaults: 64 / 10G / 03:00:00.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
''',encoding='utf-8',newline='\n')
    (P/'VERSION.txt').write_text(P.name+'\n',encoding='ascii',newline='\n')
    print(json.dumps(dict(runtime_files=len(runtime),cells=8,seeds_per_cell=1000,estimation_mode='none',
        model_weighting='maic',cpus=64,memory='10G',walltime='03:00:00',scientific_sources_modified=False)))

def package():
    files=sorted(p for p in P.rglob('*') if p.is_file() and p.relative_to(P).parts[0]
        not in ('audit','validation','__pycache__') and p.name!='SHA256SUMS')
    (P/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(P).as_posix()}\n' for p in files),encoding='ascii',newline='\n')
    archive=P.parent/(P.name+'.tar.gz')
    with tarfile.open(archive,'w:gz',compresslevel=1)as t:
        for p in files+[P/'SHA256SUMS']:t.add(p,arcname=f'{P.name}/{p.relative_to(P).as_posix()}',recursive=False)
    Path(str(archive)+'.sha256').write_text(f'{sha(archive)}  {archive.name}\n',encoding='ascii',newline='\n')
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,
        cells=8,seeds_per_cell=1000,uploaded=False,submitted=False),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',action='store_true');p.add_argument('--package',action='store_true');a=p.parse_args()
    if a.stage:stage()
    if a.package:package()
