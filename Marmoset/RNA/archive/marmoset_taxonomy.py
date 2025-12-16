import scanpy as sc
import pandas as pd
import numpy as np
import scvi

##
adata = sc.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Marmoset/HMBA_BG_Marmoset_mapped.h5ad")
adata.X = adata.raw.X.copy()

## calculate the hvg
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False, flavor="seurat_v3", batch_key="donor_name")

##
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

## Use HVGs for determining latent representation
adata_hvg = adata[:,adata.var_names[adata.var.highly_variable]].copy()
print(adata_hvg.shape)

adata_hvg.layers["UMIs"] = adata_hvg.raw[:,adata.var_names[adata.var.highly_variable]].X.copy()
 
## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    categorical_covariate_keys=["donor_name"],
)
model = scvi.model.SCVI(adata_hvg, 
                        n_hidden=128, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata.obsm["X_umap"] = sc.tl.umap(adata, min_dist=0.3, copy=True)
adata.obsm["X_umap"] = adata.obsm["X_umap"].obsm["X_umap"]

##
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/marmoset/HMBA_BG_Marmoset_mapped.h5ad")


## Assign keyword arguments (kwargs). Any unassigned args will be set to their respective function defaults.
def setup_transcriptomic_clustering(): 
    means_vars_kwargs = {
        'low_thresh': 1,
        'min_cells': 4
    }
    highly_variable_kwargs = {
        'max_genes': 3000
    }
    pca_kwargs = {
        'cell_select': 500000,
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
        # 'known_modes': known_modes_df,
        # 'similarity_threshold': 0.7
    }
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    latent_kwargs = {
        'latent_component': "X_scVI"
    }
    cluster_louvain_kwargs = {
        'k': 150, # number of nn
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 8, # cpus
        'resolution': 1.0, # resolution of louvain for taynaud method
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 50,
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.05,
            'lfc_thresh': 1.0,
            'score_thresh': 200,
            'low_thresh': 1,
            'min_genes': 5
        },
        'k': 2, # number of nn for de merge 
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
    min_samples=10,
    onestep_kwargs=onestep_kwargs,
    random_seed=123,
    tmp_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/tmp"
)

##
cluster_dict = build_cluster_dict(clusters)

##
adata.obs["cluster"] = ""
for cluster in cluster_dict.keys():
    adata.obs.cluster[cluster_dict[cluster]] = cluster

##
adata.obs.cluster = adata.obs.cluster.astype('category')

##
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalNuclei/basal_ganglia_HMBA_integrated_clustered.h5ad")
