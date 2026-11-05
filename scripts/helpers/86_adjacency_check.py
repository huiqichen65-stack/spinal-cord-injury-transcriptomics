import csv, gzip, math, os, sys
adj = sys.argv[1]
out = sys.argv[2]
opener = gzip.open if adj.endswith('.gz') else open
n = 0
tfs = set()
targets = set()
bad = 0
neg = 0
min_imp = None
max_imp = None
with opener(adj, 'rt', encoding='utf-8', newline='') as f:
    r = csv.DictReader(f, delimiter='\t')
    need = {'TF','target','importance'}
    missing = need - set(r.fieldnames or [])
    if missing:
        raise RuntimeError('Missing columns: ' + ', '.join(sorted(missing)))
    for row in r:
        n += 1
        tfs.add(row['TF'])
        targets.add(row['target'])
        try:
            imp = float(row['importance'])
        except Exception:
            bad += 1
            continue
        if not math.isfinite(imp):
            bad += 1
            continue
        if imp < 0:
            neg += 1
        if min_imp is None or imp < min_imp:
            min_imp = imp
        if max_imp is None or imp > max_imp:
            max_imp = imp
rows = [
 ('n_edges',n),
 ('n_unique_TFs',len(tfs)),
 ('n_unique_targets',len(targets)),
 ('min_importance',min_imp),
 ('max_importance',max_imp),
 ('nonfinite_importance',bad),
 ('negative_importance',neg),
 ('required_columns_present',True),
]
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['metric','value'])
    w.writerows(rows)
if n <= 0 or len(tfs) <= 0 or len(targets) <= 0 or bad > 0:
    raise RuntimeError('Adjacency QC failed')
print(f'edges={n}')
print(f'unique_TFs={len(tfs)}')
print(f'unique_targets={len(targets)}')
