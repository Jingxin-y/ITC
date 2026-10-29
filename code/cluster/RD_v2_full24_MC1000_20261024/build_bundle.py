"""Pin the final RD runtime; prepare 24 full simulation cells without running them."""
from pathlib import Path
import argparse, csv, hashlib, json, shutil, tarfile
P = Path(__file__).resolve().parent
BASE = P.parent / 'RD_final_B600_F1000_20261024' / 'code'
RR = P.parent / 'RR_v2_full24_MC1000_20261024'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p, rows):
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def stage():
    assert not (P/'VERSION.txt').exists(), 'Already staged; use --package.'
    rows = []
    for p in sorted(BASE.rglob('*')):
        if not p.is_file(): continue
        rel = p.relative_to(BASE); q = P/'code'/rel
        q.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(p, q)
        rows.append(dict(path=rel.as_posix(), md5=hashlib.md5(q.read_bytes()).hexdigest(), sha256=sha(q)))
    table(P/'RUNTIME_FILES.csv', rows)
    for name in ('env.sh', 'submit.sh', 'worker.sh', 'collect.sh', 'prepare.R'):
        s = (RR/name).read_text('utf-8')
        for a,b in [('RRFULL','RDFULL'),('RR_v2_full24','RD_v2_full24'),('RRv2','RDv2'),
                    ('RRfull','RDfull'),('formal RR','formal RD'),('Formal RR','Formal RD')]:
            s = s.replace(a,b)
        if name=='prepare.R': s=s.replace('^seed_[0-9]{6}', '^rd_seed_[0-9]{8}')
        assert 'ARDAY' not in s
        (P/name).write_text(s, encoding='utf-8', newline='\n')
    for name in ('tasks.csv','scenarios.csv','resources.csv','resource_evidence_B.csv'):
        shutil.copyfile(RR/name,P/name)
    # User resource revision 2026-10-25: only RD none cells change to 32 CPUs / 8h.
    resources=list(csv.DictReader((P/'resources.csv').open(encoding='utf-8')))
    for row in resources:
        if row['analysis']=='none':row.update(cpus=32,walltime='08:00:00')
    table(P/'resources.csv',resources)
    (P/'VERSION.txt').write_text(P.name+'\n',encoding='ascii',newline='\n')
    print(json.dumps(dict(runtime_files=len(rows),cells=24,seeds_per_cell=1000,scientific_sources_modified=False)))

def package():
    files=sorted(p for p in P.rglob('*') if p.is_file() and p.relative_to(P).parts[0]
                 not in ('audit','validation','__pycache__') and p.name!='SHA256SUMS')
    (P/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(P).as_posix()}\n' for p in files),encoding='ascii',newline='\n')
    archive=P.parent/(P.name+'.tar.gz')
    with tarfile.open(archive,'w:gz',compresslevel=1) as t:
        for p in files+[P/'SHA256SUMS']: t.add(p,arcname=f'{P.name}/{p.relative_to(P).as_posix()}',recursive=False)
    Path(str(archive)+'.sha256').write_text(f'{sha(archive)}  {archive.name}\n',encoding='ascii',newline='\n')
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,
                         cells=24,seeds_per_cell=1000,uploaded=False,submitted=False),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',action='store_true');p.add_argument('--package',action='store_true');a=p.parse_args()
    if a.stage: stage()
    if a.package: package()
