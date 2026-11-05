import csv, sys
from importlib import metadata
out = sys.argv[1]
pkgs = ['pyscenic','ctxcore','loompy','pandas','numpy','scipy','h5py']
rows = [('python', sys.version.split()[0])]
for p in pkgs:
    rows.append((p, metadata.version(p)))
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['component','version'])
    w.writerows(rows)
for k,v in rows:
    print(f'{k}={v}')
