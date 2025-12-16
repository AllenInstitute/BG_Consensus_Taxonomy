import os
import pandas as pd
import re
import sys
import matplotlib.pyplot as plt
import scipy
import numpy as np
import anndata
import h5py
import tqdm
from os.path import join as pjoin
import snapatac2 as snap
from pathlib import Path

out_path = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/ATAC'
ocs_dir = '/allen/programs/celltypes/workgroups/rnaseqanalysis/bicore/OCS_multiome_datasets/macaque/'

with h5py.File("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia/Macaque_basalganglia_AIBS_BICAN.h5ad") as f:
    anno_table = anndata.experimental.read_elem(f['obs'])

with h5py.File("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/Macaque_HMBA_basalganglia_AIT_pre-print.h5ad") as f:
    meta_table = anndata.experimental.read_elem(f['obs'])

anno_table['ar_directory'] = anno_table['ar_dir']
anno_table['load_name'] = anno_table["barcoded_cell_sample_local_name"].copy()
anno_table.loc[:,[ 'ar_directory', 'ar_id','organism']] = anno_table.loc[:,[ 'ar_directory', 'ar_id','organism']].astype(str)
anno_table['fragment_file']=anno_table.loc[:,[ 'ar_directory', 'ar_id']].apply(lambda x: '/'.join(x),axis=1)+'/outs/atac_fragments.tsv.gz'

def find_atac_fragments(base_dir):
    base_path = Path(base_dir)
    records = []
    for file in base_path.rglob("atac_fragments.tsv.gz"):
        parent_parts = file.parent.parts
        folder_minus2 = parent_parts[-2] if len(parent_parts) >= 2 else None
        folder_minus3 = parent_parts[-3] if len(parent_parts) >= 3 else None
        records.append({
            "file_path": str(file),
            "mystery_hash": folder_minus2,
            "load_name": folder_minus3
        })
    return pd.DataFrame(records)

df = find_atac_fragments(ocs_dir)
df_frag_dict = dict(zip(df['load_name'], df['file_path']))
anno_table.loc[anno_table['load_name'].isin(df['load_name']), 'fragment_file'] = anno_table.loc[anno_table['load_name'].isin(df['load_name']),'load_name'].replace(df_frag_dict)

atac_samples = anno_table.loc[:, ['organism', 'fragment_file', 'ar_directory', 'ar_id']].drop_duplicates().reset_index()
atac_samples = atac_samples.loc[~atac_samples['fragment_file'].str.contains('nan/nan'), :]

## Write out
atac_samples.to_csv(pjoin(out_path, 'atac_samples.csv'), index=False)

## Merge atac_samples and meta_table
missing_columns = list(set(anno_table.columns) - set(meta_table.columns))
anno = pd.merge(anno_table.loc[:,missing_columns], meta_table, how='inner', left_index=True, right_index=True)
anno.to_csv(pjoin(out_path, 'anno_table.csv'))
