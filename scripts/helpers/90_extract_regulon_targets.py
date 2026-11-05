import csv
import sys
from pyscenic.cli.utils import load_signatures

reg_path = sys.argv[1]
tf_path = sys.argv[2]
out_path = sys.argv[3]

with open(tf_path, 'rt', encoding='utf-8') as f:
    tfs = [x.strip() for x in f if x.strip()]

wanted = set(tfs)
sigs = list(load_signatures(reg_path))

rows = []
found = set()

for s in sigs:
    reg_name = str(s.name)
    tf = reg_name
    if tf.endswith('(+)') or tf.endswith('(-)'):
        tf = tf[:-3]
    tf = tf.strip()

    if tf not in wanted:
        continue

    found.add(tf)
    gene2weight = dict(s.gene2weight)
    genes = list(gene2weight.keys())

    for g in genes:
        wgt = gene2weight[g]
        rows.append((tf, reg_name, str(g), float(wgt), len(genes)))

missing = [tf for tf in tfs if tf not in found]
if missing:
    raise RuntimeError('Top TF regulons not found: ' + ', '.join(missing))

with open(out_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['TF','regulon','target','regulon_weight','regulon_size'])
    w.writerows(rows)

print(f'top_tfs={len(found)}')
print(f'target_rows={len(rows)}')
