import csv
import os
import sys
import numpy as np
import pandas as pd

auc_path = sys.argv[1]
cell_path = sys.argv[2]
final_gz = sys.argv[3]
qc_path = sys.argv[4]

import gzip

def read_lines(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        return [x.rstrip('\r\n') for x in f if x.rstrip('\r\n') != '']

expected_cells = read_lines(cell_path)

df = pd.read_csv(auc_path, index_col=0)

df.index = df.index.astype(str)
values = df.to_numpy(dtype=float, copy=False)

cell_order_match = df.index.tolist() == expected_cells
cell_set_match = set(df.index.tolist()) == set(expected_cells)
unique_cells = df.index.is_unique
unique_regulons = df.columns.is_unique

finite = np.isfinite(values)
nonfinite = int(values.size - finite.sum())
negative = int((values < 0).sum())

col_min = np.nanmin(values, axis=0)
col_max = np.nanmax(values, axis=0)
zero_variance = int(np.sum(np.isclose(col_max, col_min)))
all_zero = int(np.sum(np.all(np.isclose(values, 0.0), axis=0)))

rows = [
    ('n_cells', df.shape[0]),
    ('n_regulons', df.shape[1]),
    ('cell_order_match', cell_order_match),
    ('cell_set_match', cell_set_match),
    ('unique_cells', unique_cells),
    ('unique_regulons', unique_regulons),
    ('nonfinite_values', nonfinite),
    ('negative_values', negative),
    ('zero_variance_regulons', zero_variance),
    ('all_zero_regulons', all_zero),
    ('global_min_AUC', float(np.nanmin(values))),
    ('global_max_AUC', float(np.nanmax(values))),
]

with open(qc_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['metric','value'])
    w.writerows(rows)

if not cell_order_match:
    raise RuntimeError('AUCell cell order does not match the full atlas.')

if nonfinite > 0 or negative > 0:
    raise RuntimeError('AUCell matrix contains invalid values.')

if all_zero > 0:
    raise RuntimeError('One or more regulons are all-zero.')

df.to_csv(final_gz, compression='gzip')

print(f'AUCell QC PASS: {df.shape[0]} cells x {df.shape[1]} regulons')
