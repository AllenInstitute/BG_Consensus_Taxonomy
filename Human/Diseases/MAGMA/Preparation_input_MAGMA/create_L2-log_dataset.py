'''
Code is adapted from Duncan et al. (2024)
Create the level 2 (cluster level) Siletti matrix. Cellular expression value is ln(1+p) transformed prior to aggregation.
'''
import h5py
import os
import numpy as np
import numexpr as ne # to speed up
from tqdm import tqdm
import anndata
import scipy.sparse
import pandas as pd

level = 'Cluster'
level_name = 'cluster_id'

base_dir_BG_dataset = ".../data/newest_BG/"
new_file_name = paste0(base_dir_BG_dataset, "newest_BG_log2.h5")

adata = anndata.read_h5ad('/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/Human_basalganglia_HMBA_AIT19-5_anno_latest.h5ad')

df_cluster = adata.obs[level_name].to_numpy()

n_genes = adata.shape[1]
n_clusters = len(np.unique(df_cluster))


f_write = h5py.File(new_file_name, 'w')
avg_matrix = f_write.create_dataset("matrix", (n_genes, n_clusters), '<f8')

for cluster_num in tqdm(np.arange(n_clusters)):
    cluster_str = str(cluster_num)
    dset_sub = adata.raw[adata.obs[level_name] == cluster_str, :]

    if scipy.sparse.issparse(dset_sub.X):
        dset_sub = dset_sub.X.toarray()  # Convert sparse to dense
    else:
        dset_sub = dset_sub.X

    dset_sub = dset_sub.astype(np.float64)
    test_size = ne.evaluate('sum(log1p(dset_sub), axis=0)')/dset_sub.shape[0]
    print(test_size.shape)
    avg_matrix[:,cluster_num] = test_size

# Create labels
c_f_write = f_write.create_dataset(level, data=clusters)

dt = h5py.string_dtype(encoding='utf-8')
gene_names = np.array(adata.var["gene_id"].values, dtype=dt)
f_write.create_dataset("gene_ids", data=gene_names)

f_write.close()
