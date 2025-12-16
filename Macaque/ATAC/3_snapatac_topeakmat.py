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
import pybedtools
import anndata as ad

root_directory = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/ATAC'
groupby = 'Group'

##
file_path = os.path.join(root_directory, 'h5ad/concatenated_102225.h5ads')
data = snap.read_dataset(file_path, mode='r+')
data.obs[groupby] = data.adatas.obs[groupby]

##
peak_path = os.path.join(root_directory, 'merged_peaks.bed')

##
peak_adata = snap.pp.make_peak_matrix(data, peak_file=peak_path)

## Write out
peak_adata.write_h5ad(os.path.join(root_directory, 'h5ad/concatenated_peakmat.h5ad'))

##
data.close()

##
snap.pp.select_features(peak_adata, n_features=250000)
snap.tl.spectral(peak_adata)
snap.tl.umap(peak_adata)

##
snap.pl.umap(peak_adata, color='Group', interactive=False, height=500)

## Add donor
meta_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
donor_meta = pd.read_csv(os.path.join(meta_dir, 'consensus_hmba_basalganglia_metadata.csv'), index_col=0)

human_donor_meta = donor_meta[donor_meta['organism']=='Human']
donor_dict = dict(zip(human_donor_meta['load_id'], human_donor_meta['donor_id']))

## Gather sample info from cell barcodes which are barcode-load_id
peak_adata.obs["load_id"] = [
    re.sub(r'^.*?-', '', cell_id) for cell_id in peak_adata.obs_names
]
peak_adata.obs["donor"] = peak_adata.obs["load_id"].map(donor_dict)
peak_adata.write_h5ad(os.path.join(root_directory, 'h5ad/concatenated_peakmat_UMAP.h5ad'))

## Harmony batch correction
peak_adata = ad.read_h5ad(os.path.join(root_directory, 'h5ad/concatenated_peakmat_UMAP.h5ad'))
snap.pp.harmony(peak_adata, batch="donor", max_iter_harmony=20)

snap.tl.umap(peak_adata, use_rep="X_spectral_harmony")

umap_df = pd.DataFrame(peak_adata.obsm['X_umap'], index=peak_adata.obs_names, columns=['UMAP1', 'UMAP2'])
umap_df["Group"] = peak_adata.obs["Group"].values
umap_df.to_csv(os.path.join(root_directory, 'h5ad/concatenated_peakmat_UMAP_harmony.csv'))

## Write out
peak_adata.write_h5ad(os.path.join(root_directory, 'h5ad/concatenated_peakmat_UMAP.h5ad'))
peak_adata.close()