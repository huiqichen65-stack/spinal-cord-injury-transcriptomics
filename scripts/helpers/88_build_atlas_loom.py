import csv
import gzip
import os
import sys
import numpy as np
import loompy
from scipy.io import mmread

mtx_path = sys.argv[1]
gene_path = sys.argv[2]
cell_path = sys.argv[3]
loom_path = sys.argv[4]
qc_path = sys.argv[5]

def read_lines(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        return [x.rstrip('\r\n') for x in f if x.rstrip('\r\n') != '']

genes = read_lines(gene_path)
cells = read_lines(cell_path)

mat = mmread(mtx_path).tocsc()

if mat.shape != (len(genes), len(cells)):
    raise RuntimeError(
        f'Matrix shape {mat.shape} does not match genes/cells '
        f'({len(genes)}, {len(cells)})'
    )

if os.path.exists(loom_path):
    os.remove(loom_path)

chunk = 500

row_attrs = {
    'Gene': np.asarray(genes, dtype=object)
}

with loompy.new(loom_path) as ds:
    for start in range(0, len(cells), chunk):
        end = min(start + chunk, len(cells))
        block = mat[:, start:end].toarray().astype(np.float32, copy=False)
        col_attrs = {
            'CellID': np.asarray(cells[start:end], dtype=object)
        }
        if start == 0:
            ds.add_columns(
                block,
                col_attrs=col_attrs,
                row_attrs=row_attrs
            )
        else:
            ds.add_columns(
                block,
                col_attrs=col_attrs
            )
        print(f'loom columns: {end}/{len(cells)}', flush=True)

with loompy.connect(loom_path, mode='r') as ds:
    loom_genes = [str(x) for x in ds.ra.Gene]
    loom_cells = [str(x) for x in ds.ca.CellID]
    shape = ds.shape

gene_order_match = loom_genes == genes
cell_order_match = loom_cells == cells
shape_match = shape == (len(genes), len(cells))

rows = [
    ('n_genes', shape[0]),
    ('n_cells', shape[1]),
    ('shape_match', shape_match),
    ('gene_order_match', gene_order_match),
    ('cell_order_match', cell_order_match),
    ('loom_size_bytes', os.path.getsize(loom_path)),
]

with open(qc_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['metric','value'])
    w.writerows(rows)

if not (shape_match and gene_order_match and cell_order_match):
    raise RuntimeError('Full-atlas loom QC failed.')

print('Full-atlas loom QC: PASS', flush=True)
