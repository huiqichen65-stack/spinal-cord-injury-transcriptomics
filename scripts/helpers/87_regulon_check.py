import csv, os, sys
from pyscenic.cli.utils import load_signatures

reg_file = sys.argv[1]
qc_out = sys.argv[2]
summary_out = sys.argv[3]

sigs = list(load_signatures(reg_file))

if len(sigs) == 0:
    raise RuntimeError('No regulons/signatures could be loaded from ctx output.')

rows = []
all_targets = set()

for s in sigs:
    name = str(s.name)
    genes = list(s.genes)
    all_targets.update(genes)
    rows.append((name, len(genes)))

names = [x[0] for x in rows]
n_positive = sum('(+' in n or '(+)' in n for n in names)
n_negative = sum('(-' in n or '(-)' in n for n in names)

with open(summary_out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['regulon','n_genes'])
    w.writerows(rows)

qc = [
    ('regulon_file_exists', os.path.exists(reg_file)),
    ('regulon_file_size_bytes', os.path.getsize(reg_file)),
    ('n_regulons', len(sigs)),
    ('n_unique_regulon_names', len(set(names))),
    ('n_unique_target_genes', len(all_targets)),
    ('min_regulon_size', min(x[1] for x in rows)),
    ('median_regulon_size', sorted(x[1] for x in rows)[len(rows)//2]),
    ('max_regulon_size', max(x[1] for x in rows)),
    ('n_positive_named', n_positive),
    ('n_negative_named', n_negative),
]

with open(qc_out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['metric','value'])
    w.writerows(qc)

print(f'n_regulons={len(sigs)}')
print(f'n_unique_targets={len(all_targets)}')
