import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import pickle
import scvi
import os

import transcriptomic_clustering as tc
from transcriptomic_clustering.iterative_clustering import (
    build_cluster_dict, iter_clust, OnestepKwargs, summarize_final_clusters, 
)
from transcriptomic_clustering.final_merging import (
    final_merge, FinalMergeKwargs,
)

##
species = "Macaque"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
data_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
figure_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia/figures"
model_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia/models"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

os.chdir(work_dir)

## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT11-8_anno_filtered_scVI.h5ad"))

##
def setup_transcriptomic_clustering():
    means_vars_kwargs = {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, 0.6931472 to match to bigcat
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    }
    highly_variable_kwargs = {
        'max_genes': 4000 # originally 3000, 4000 to match to bigcat
    }
    pca_kwargs = {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    }
    ## Leave empty if you don't want to use known_modes
    filter_known_modes_kwargs = {
        # 'known_modes': known_modes_df, # a pd dataframe. index is obs (cell) names, columns are known modes. Originally commented out
        'similarity_threshold': 0.7
    }
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    latent_kwargs = {
        'latent_component': "X_scVI"
    }
    cluster_louvain_kwargs = {
        'k': 15, # number of nn, originally 150, change to 15
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus # originally 8`
        'resolution': 1.0 # resolution of louvain for taynaud method
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, 10 to match to bigcat
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.05,
            'lfc_thresh': 0.6931472, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, 100 to match to bigcat
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, 0.6931472 to match to bigcat
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, 4 to match to bigcat
        'de_method': 'ebayes'
    }
 
    onestep_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )
    return onestep_kwargs

##
onestep_kwargs = setup_transcriptomic_clustering()

## Run clustering
clusters, markers = iter_clust(
    adata,
    min_samples=4,
    onestep_kwargs=onestep_kwargs,
    random_seed=123,
    tmp_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/analysis/tmp"
)

##
with open(os.path.join(work_dir, 'clustering', 'clustering_results.pkl'), 'wb') as f:
    pickle.dump(clusters, f)

##
with open(os.path.join(work_dir, 'clustering', 'clustering_markers.pkl'), 'wb') as f:
    pickle.dump(markers, f)

## -------------------------
## Final merging using all clusters

## Convert to a list of lists
clusters_as_lists = [array.tolist() for array in clusters]

def setup_merging(): 
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, 
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 0.6931472, 
            'score_thresh': 100, 
            'low_thresh': 0.6931472, 
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
    }
    latent_kwargs = {  
        "latent_component": "X_scVI"
    }
    merge_kwargs = FinalMergeKwargs(
        merge_clusters_kwargs=merge_clusters_kwargs,
        latent_kwargs=latent_kwargs  
    )
    return merge_kwargs

merge_kwargs = setup_merging()

## Run the final merging
clusters_after_merging, markers = final_merge(
    adata, 
    clusters_as_lists, 
    # markers, # required for PCA, but optional if using a pre-computed latent space
    n_samples_per_clust=20, 
    random_seed=123, 
    n_jobs = 30, # modify this to the number of cores you want to use
    return_markers_df = False, # return the pair-wise DE results for each cluster pair. If False (default), only return a set of markers (top 20 of up and down regulated genes in each pair comparison)
    final_merge_kwargs=merge_kwargs
)

## Gather clustering resutls
n_cells = sum(len(i) for i in clusters_after_merging)
cl = ['unknown']*n_cells
for i in range(len(clusters_after_merging)):
    for j in clusters_after_merging[i]:
        cl[j] = i+1

##
clusters = pd.DataFrame({'cl': cl}, index=adata.obs_names)
np.all(adata.obs_names == clusters.index)

## Save the final clustering results
adata.obs['cluster'] = clusters['cl']
adata.obs['cluster'] = adata.obs['cluster'].astype("category")
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9.h5ad")) 
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_HMBA_AIT11-9.h5ad")) 