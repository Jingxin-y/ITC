from pathlib import Path
import csv,hashlib,shutil,tarfile,json
P=Path(__file__).resolve().parent
BASE=P.parent/'RR_effect_v2_20261015'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def stage():
    assert not (P/'VERSION.txt').exists(),'Already staged; use --package.'
    rows=[]
    for folder in ('data','estimation','summary'):
        for p in sorted((BASE/folder).rglob('*')):
            if not p.is_file() or p.suffix.lower() not in ('.r','.md','.txt','.csv','.json'):continue
            if any(x in ('validation','output','.backups') for x in p.relative_to(BASE).parts):continue
            rel=p.relative_to(BASE);q=P/'code'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
            rows.append(dict(path=rel.as_posix(),md5=hashlib.md5(q.read_bytes()).hexdigest(),sha256=sha(q)))
    for name in ('load_anchored.R','BRMPLUS_LICENSE.txt'):
        q=P/'code'/name;shutil.copyfile(BASE/name,q)
        rows.append(dict(path=name,md5=hashlib.md5(q.read_bytes()).hexdigest(),sha256=sha(q)))
    assert (P/'code/estimation/helpers/CI_marginal_target.R').read_bytes()==(P.parent/'RR_final_B600_F1000_20261024/code/estimation/helpers/CI_marginal_target.R').read_bytes()
    assert 'webster2024_rr_ess_em_v2' in (P/'code/data/generate_webster_RR_scenarios.R').read_text('utf-8')
    table(P/'RUNTIME_FILES.csv',rows)
    cases=[]
    for event in ('common','rare'):
        for n in (100,500):
            for ess in (.15,.5):cases.append(dict(scenario_id=len(cases)+1,event=event,n=n,ess_ratio=ess))
    table(P/'scenarios.csv',cases)
    tasks=[]
    for c in cases:
        for analysis,est,weight in (('none','none','unit'),('score_unit','scores','unit'),('score_MAIC','scores','maic')):
            tasks.append(dict(task_id=len(tasks)+1,**c,analysis=analysis,estimation_mode=est,model_weighting=weight,
                ess_source='diffuse',em_level='high',weight_mode='all',nsim=1000,first_seed=1,
                cell=f"S{c['scenario_id']:02d}_{c['event']}_n{c['n']}_ess{round(100*c['ess_ratio']):03d}_{analysis}"))
    table(P/'tasks.csv',tasks)
    table(P/'resources.csv',[dict(analysis=analysis,n=n,cpus=16 if analysis=='none'else 96,
        memory='4G' if analysis=='none'else '20G',walltime='04:00:00'if analysis=='none'else
        ('18:00:00'if n==100 else '19:00:00'))for analysis in ('none','score_unit','score_MAIC')for n in (100,500)])
    (P/'VERSION.txt').write_text(P.name+'\n',encoding='ascii',newline='\n')
    print('Staged v2 generator + final RR CI, 24 cells x 1000 seeds; no scientific source edited.')
def package():
    files=sorted(p for p in P.rglob('*') if p.is_file() and p.relative_to(P).parts[0] not in ('audit','validation','__pycache__') and p.name!='SHA256SUMS')
    (P/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(P).as_posix()}\n' for p in files),encoding='ascii',newline='\n')
    arc=P.parent/(P.name+'.tar.gz')
    with tarfile.open(arc,'w:gz',compresslevel=1) as t:
        for p in files+[P/'SHA256SUMS']:t.add(p,arcname=f'{P.name}/{p.relative_to(P).as_posix()}',recursive=False)
    Path(str(arc)+'.sha256').write_text(f'{sha(arc)}  {arc.name}\n',encoding='ascii',newline='\n')
    print(json.dumps(dict(archive=str(arc),sha256=sha(arc),bytes=arc.stat().st_size,compute_tasks=24,seeds_per_task=1000,submitted=False),indent=2))
if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--stage',action='store_true');a.add_argument('--package',action='store_true');v=a.parse_args()
    if v.stage:stage()
    if v.package:package()
