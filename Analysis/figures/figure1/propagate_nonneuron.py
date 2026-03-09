import os, sys
import seaborn
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import tqdm
import numpy as np

##
adata = sc.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_non-neurons.h5ad")

from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
from sklearn.neighbors import KNeighborsClassifier

def distribute_evenly(a, k):
    n = len(a)  # Number of elements
    r = [0] * n  # Initialize result vector with zeros
    # First pass: Try to distribute evenly based on the initial total
    avg = k // n
    for i in range(n):
        r[i] = min(a[i], avg)
    remaining = k - sum(r)  # Calculate remaining after initial distribution
    # Distribute the remaining as evenly as possible
    while remaining > 0:
        # Calculate the updated average for the remaining amount
        updates = 0
        for i in range(n):
            if r[i] < a[i]:
                potential_add = min(a[i] - r[i], remaining // (n - updates))
                if potential_add > 0:
                    r[i] += potential_add
                    remaining -= potential_add
                    updates += 1
                    if remaining == 0:
                        break
        # In case there are no updates, but still remaining (to avoid infinite loop)
        if updates == 0:
            break
    return r


def select_k_cells(adata, obs_column, k):
    a = adata.obs[obs_column].value_counts()
    categories = a.index
    r = distribute_evenly(a,k)
    inds=[]
    for i,cat in enumerate(categories):
        cat_cells = adata.obs.index[adata.obs[obs_column] == cat].tolist()
        inds.append(np.random.choice(cat_cells, size=r[i], replace=False))
    selected_indices=[item for sublist in inds for item in sublist]
    return adata[selected_indices,:]

def consensus_knn_transfer(adata, label_col, train_mask, test_mask, n_neighbors=20, n_iter=10):
    preds_list = []
    adata.obs[label_col + '_propagated'] = list(adata.obs[label_col])
    train_data = adata[train_mask]
    # Determine number of cells to sample from each class
    min_count = train_data.obs[label_col].value_counts().min()
    k_val = len(train_data.obs[label_col].unique()) * min_count * 2
    for _ in tqdm.tqdm(range(n_iter)):
        # Downsample training data evenly
        equalized_train = select_k_cells(train_data, label_col, k_val)
        neigh = KNeighborsClassifier(n_neighbors=n_neighbors)
        X_train = equalized_train.obsm['X_scVI']
        y_train = equalized_train.obs[label_col]
        neigh.fit(X_train, y_train)
        test_data = adata[test_mask]
        X_test = test_data.obsm['X_scVI']
        preds = neigh.predict(X_test)
        preds_list.append(preds)
    # Stack predictions from all iterations and take majority vote for each test cell.
    preds_array = np.vstack(preds_list)  # shape: (n_iter, n_test)
    consensus_preds = [
        Counter(preds_array[:, j]).most_common(1)[0][0]
        for j in range(preds_array.shape[1])
    ]
    adata.obs.loc[test_mask, label_col + '_propagated'] = consensus_preds
    return adata


adata.obs['Group'] = adata.obs['Group'].astype(str)

train_mask = adata.obs['Group']!='nan'
test_mask  = adata.obs['Group']=='nan'
adata = consensus_knn_transfer(adata, 'Group', train_mask, test_mask,
                               n_neighbors=30, n_iter=20)

##
adata.write('/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_non-neurons_labeltransfer.h5ad')

#####################################
## Update with all genes
#####################################
adata = sc.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v7.h5ad")

import h5py
import anndata
with h5py.File("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_non-neurons_labeltransfer.h5ad") as f:
    obsm = anndata.experimental.read_elem(f['obsm'])
    obs = anndata.experimental.read_elem(f['obs'])
    var = anndata.experimental.read_elem(f['var'])

for k in obsm.keys():
    obsm[k] = pd.DataFrame(obsm[k],index=obs.index)

adata = adata[~adata.obs.index.duplicated(),:]
obs = obs.loc[~obs.index.duplicated()]
for k in obsm.keys():
    obsm[k] = obsm[k].loc[~obsm[k].index.duplicated()]

adata.obs.index.duplicated().sum()
obs.index.duplicated().sum()

obs = obs.loc[adata.obs.index[adata.obs.index.isin(obs.index)]]
adata = adata[[str(x) for x in obs.index], :]
obs = obs.loc[adata.obs_names]

adata.obs = adata.obs.copy()

for col in ['Group_propagated']:
    print(col)
    adata.obs[col] = list(obs[col].astype(obs[col].dtype))

for k in obsm.keys():
    adata.obsm[k] = obsm[k].loc[adata.obs_names].to_numpy()

adata.obs

adata.obs['Subclass_propagated'] = adata.obs['Group_propagated'].replace(adata.obs.groupby('Group')['Subclass'].value_counts().unstack().idxmax(1).to_dict())
adata.obs['Class_propagated'] = adata.obs['Subclass_propagated'].replace(adata.obs.groupby('Subclass')['Class'].value_counts().unstack().idxmax(1).to_dict())

adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_non-neuron_alignment_labeltransfer_allgenes.h5ad")
adata.write('/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_non-neurons_labeltransfer_allgenes.h5ad')
