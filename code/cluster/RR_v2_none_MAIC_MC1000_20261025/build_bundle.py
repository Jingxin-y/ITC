from pathlib import Path
import argparse,csv,hashlib,json,shutil,tarfile
P=Path(__file__).resolve().parent
BASE=P.parent/'RR_v2_full24_MC1000_20261024'
AUTO=P.parent/'RR_none_autoresume_20261025'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,s):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s,encoding='utf-8',newline='\n')
def table(p,rows):
    with p.open('w',encoding='utf-8',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def stage():
    assert not(P/'VERSION.txt').exists(),'Already staged; use --package.'
    shutil.copytree(BASE/'code',P/'code')
    shutil.copyfile(BASE/'RUNTIME_FILES.csv',P/'RUNTIME_FILES.csv')
    shutil.copyfile(BASE/'scenarios.csv',P/'scenarios.csv')
    old=list(csv.DictReader((BASE/'tasks.csv').open(encoding='utf-8')))
    rows=[]
    for r in old:
        if r['analysis']!='none':continue
        r['task_id']=len(rows)+1;r['analysis']='none_MAIC';r['model_weighting']='maic'
        r['cell']=r['cell']+'_MAIC';rows.append(r)
    assert len(rows)==8
    table(P/'tasks.csv',rows)
    table(P/'resources.csv',[dict(analysis='none_MAIC',n=n,cpus=64,memory='10G',walltime='04:00:00')for n in(100,500)])
    common=(BASE/'common.R').read_text('utf-8')
    i=common.index('tasks<-function(home){');j=common.index('set.analysis<-function(r){')
    common=common[:i]+'''tasks<-function(home){
    x<-read.csv(file.path(home,'tasks.csv'),stringsAsFactors=FALSE)
    must(nrow(x)==8L&&identical(x$task_id,1:8)&&identical(x$scenario_id,1:8)&&
        all(x$nsim==1000L)&&all(x$first_seed==1L)&&all(x$analysis=='none_MAIC')&&
        all(x$estimation_mode=='none')&&all(x$model_weighting=='maic')&&
        all(x$ess_source=='diffuse')&&all(x$em_level=='high')&&all(x$weight_mode=='all')&&
        identical(x$event,rep(c('common','rare'),each=4L))&&
        identical(x$n,rep(c(100L,100L,500L,500L),2L))&&
        isTRUE(all.equal(x$ess_ratio,rep(c(.15,.5),4L)))&&!anyDuplicated(x$cell),
        'Invalid eight-case none+MAIC design')
    x
}
''' +common[j:]
    i=common.index('resources<-function(home,r){');j=common.index('check.config<-function(d,cfg){')
    common=common[:i]+'''resources<-function(home,r){
    z<-read.csv(file.path(home,'resources.csv'),stringsAsFactors=FALSE)
    z<-z[z$analysis==r$analysis&z$n==r$n,,drop=FALSE]
    must(nrow(z)==1L&&z$cpus==64L&&z$memory=='10G'&&z$walltime=='04:00:00','Expected 64 CPUs / 10G / 4h')
    z
}
''' +common[j:]
    write(P/'common.R',common)
    for name in('prepare.R','run.R','collect.R','submit.sh','worker.sh','collect.sh'):
        s=(BASE/name).read_text('utf-8').replace('RRFULL_','RRNMW_')
        s=s.replace('1:24','1:8').replace('1..24','1..8').replace('ALL 24 CELLS','ALL 8 NONE+MAIC CELLS')
        s=s.replace('RRv2_1000_g','RRnoneMAIC_g').replace('RRv2_report','RRnoneMAIC_report').replace('RRfull_%A_%a','RRnoneMAIC_%A_%a')
        write(P/name,s)
    write(P/'env.sh','''#!/bin/bash
export RRNMW_HOME
RRNMW_HOME=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export RRNMW_OUT="${RRNMW_OUT:-${SCRATCH:-/scratch/$USER}/RR_v2_none_MAIC_MC1000_results_20261025}"
export RRNMW_ACCOUNT="${RRNMW_ACCOUNT:-def-liteep}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
''')
    for name in('inspect_none.R','submit_auto_none.sh','controller.sh'):
        s=(AUTO/name).read_text('utf-8').replace('RRFULL_','RRNMW_').replace('RRNONE_','RRNMW_AUTO_')
        s=s.replace('RR_v2_full24_MC1000_20261024','RR_v2_none_MAIC_MC1000_20261025')
        s=s.replace('${RRNMW_AUTO_BASE:-$HOME/RR_v2_none_MAIC_MC1000_20261025}',
            '${RRNMW_AUTO_BASE:-$(cd -- "$auto_home/.." && pwd -P)}')
        s=s.replace('RRv2_1000_g','RRnoneMAIC_g').replace('RRnone_check','RRnoneMAIC_check').replace('RRnone_resume','RRnoneMAIC_resume').replace('RRnone_finalreport','RRnoneMAIC_finalreport')
        s=s.replace('auto_none_20261025','auto_none_MAIC_20261025')
        s=s.replace("tab$analysis=='none'","tab$analysis=='none_MAIC'")
        s=s.replace('c(1L,4L,7L,10L,13L,16L,19L,22L)','1:8')
        s=s.replace("tab$model_weighting=='unit'","tab$model_weighting=='maic'")
        s=s.replace("identical(s$model_weighting,'unit')","identical(s$model_weighting,'maic')")
        s=s.replace('1|4|7|10|13|16|19|22','1|2|3|4|5|6|7|8')
        s=s.replace('--cpus-per-task=16 --mem=4G --time=04:00:00','--cpus-per-task=64 --mem=10G --time=04:00:00')
        write(P/'auto'/name,s)
    write(P/'submit_auto.sh','''#!/bin/bash
set -euo pipefail
[[ "$#" == 0 ]] || { echo 'Usage: bash submit_auto.sh (all eight none+MAIC cases)'; exit 1; }
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$here/submit.sh"
bash "$here/auto/submit_auto_none.sh"
''')
    write(P/'VERSION.txt',P.name+'\n')
    print('STAGED: eight none+MAIC cases, original scientific source unchanged, 64 CPUs /10G /4h.')
def package():
    af=sorted(x for x in(P/'auto').iterdir()if x.is_file()and x.name!='SHA256SUMS')
    write(P/'auto/SHA256SUMS',''.join(f'{sha(x)}  {x.name}\n'for x in af))
    files=sorted(x for x in P.rglob('*')if x.is_file()and x.relative_to(P).parts[0]not in('audit','validation','__pycache__')and x!=P/'SHA256SUMS')
    write(P/'SHA256SUMS',''.join(f'{sha(x)}  {x.relative_to(P).as_posix()}\n'for x in files))
    arc=P.parent/(P.name+'.tar.gz')
    with tarfile.open(arc,'w:gz',compresslevel=1)as t:
        for x in files+[P/'SHA256SUMS']:t.add(x,arcname=f'{P.name}/{x.relative_to(P).as_posix()}',recursive=False)
    write(Path(str(arc)+'.sha256'),f'{sha(arc)}  {arc.name}\n')
    print(json.dumps(dict(archive=str(arc),sha256=sha(arc),bytes=arc.stat().st_size,cases=8,seeds_per_case=1000,
        cpus=64,memory='10G',time='04:00:00',submitted=False),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',action='store_true');p.add_argument('--package',action='store_true');a=p.parse_args()
    if a.stage:stage()
    if a.package:package()
