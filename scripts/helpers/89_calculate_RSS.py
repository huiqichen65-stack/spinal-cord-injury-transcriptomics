import sys
import numpy as np
import pandas as pd
from pyscenic.rss import regulon_specificity_scores

auc_path = sys.argv[1]
metadata_path = sys.argv[2]
state_col = sys.argv[3]
rss_out = sys.argv[4]

auc = pd.read_csv(auc_path, index_col=0)
meta = pd.read_csv(metadata_path)

auc.index = auc.index.astype(str)
meta['cell'] = meta['cell'].astype(str)

if state_col not in meta.columns:
    raise RuntimeError(f'Metadata column not found: {state_col}')

if meta['cell'].duplicated().any():
    raise RuntimeError('Duplicate cell IDs in metadata.')

if auc.index.duplicated().any():
    raise RuntimeError('Duplicate cell IDs in AUCell matrix.')

if set(auc.index) != set(meta['cell']):
    raise RuntimeError('AUCell and metadata cell sets differ.')


meta = meta.set_index('cell').loc[auc.index]

if meta[state_col].isna().any():
    raise RuntimeError('NA state labels found.')

labels = meta[state_col].astype(str)

if (auc.to_numpy(dtype=float) < 0).any():
    raise RuntimeError('Negative AUC value detected.')

if not np.isfinite(auc.to_numpy(dtype=float)).all():
    raise RuntimeError('Non-finite AUC value detected.')

if (auc.sum(axis=0) <= 0).any():
    raise RuntimeError('At least one regulon has zero total AUC.')

rss = regulon_specificity_scores(auc, labels)

if not np.isfinite(rss.to_numpy(dtype=float)).all():
    raise RuntimeError('Non-finite RSS value detected.')

if (rss.to_numpy(dtype=float) < 0).any() or (rss.to_numpy(dtype=float) > 1).any():
    raise RuntimeError('RSS value outside [0,1].')

rss.to_csv(rss_out)

print(f'cells={auc.shape[0]}')
print(f'regulons={auc.shape[1]}')
print(f'states={rss.shape[0]}')
print(f'rss_min={rss.to_numpy().min()}')
print(f'rss_max={rss.to_numpy().max()}')
