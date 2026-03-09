from __future__ import annotations

import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
import h5py
import os
import csv
import re
from scipy.sparse import csr_matrix

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
out_dir = os.path.join(analysis_dir, "atac", "specificity")

## -------------------------------------------
## Load Macaque chromosome alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## AnnoTable
cluster_meta = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## From CRESted
species_h5ads = {
    "human" : os.path.join(out_dir, "human_with_gini.h5ad"),
    "macaque" : os.path.join(out_dir, "macaque_with_gini.h5ad"),
    "marmoset": os.path.join(out_dir, "marmoset_with_gini.h5ad"),
}

for species, h5ad_path in species_h5ads.items():
    print(f"Processing {species}...")
    species_adata = ad.read_h5ad(h5ad_path)
    del species_adata.obs["file_path"]
    ## Identify most accessible cell type for each peak
    max_idx = np.argmax(species_adata.X, axis=0)
    max_cell_types = species_adata.obs_names[max_idx]
    species_adata.var["max_cell_type"] = max_cell_types.tolist()
    ## Gather top 100 peaks per cell type based on specificity (Gini)
    top_peaks = []
    for cell_type in species_adata.obs_names.unique():
        cell_type_peaks = species_adata.var[species_adata.var["max_cell_type"] == cell_type]
        top_cell_type_peaks = cell_type_peaks.nlargest(100, "gini_scores")
        top_peaks.append(top_cell_type_peaks.index.tolist())
    ##
    top_peaks = [peak for sublist in top_peaks for peak in sublist]  # Flatten list
    top_peaks = list(set(top_peaks))  # Unique peaks
    ## Subset for heatmap
    heatmap_adata = species_adata[:,top_peaks].copy()
    heatmap_adata.obs["cell_type"] = heatmap_adata.obs_names.astype("category")
    ## Clear non-essential layers/obsm/uns
    heatmap_adata.layers.clear()
    heatmap_adata.obsm.clear()
    heatmap_adata.uns.clear()
    ## Set nan to 0
    heatmap_adata.X = np.nan_to_num(heatmap_adata.X, nan=0.0)
    ## Normalize for visualization
    sc.pp.normalize_total(heatmap_adata, target_sum=1e4)
    sc.pp.log1p(heatmap_adata)
    # sc.pp.scale(heatmap_adata)
    ## Arrange by max_cell_type var and Group order
    group_order = cluster_meta.Group
    ## If any Groups are missing add a obs of nans
    missing_groups = [grp for grp in group_order if grp not in heatmap_adata.obs_names]
    for missing in missing_groups:
        missing_adata = ad.AnnData(
            X=np.full((1, heatmap_adata.n_vars), np.nan),
            obs=pd.DataFrame(index=[missing]),
            var=heatmap_adata.var.copy()
        )
        missing_adata.obs["cell_type"] = missing
        heatmap_adata = ad.concat([heatmap_adata, missing_adata], axis=0, join="outer", merge="same")
    ## Ensure correct categorical ordering for obs
    heatmap_adata = heatmap_adata[heatmap_adata.obs_names.isin(group_order),:]
    heatmap_adata.obs["cell_type"] = heatmap_adata.obs["cell_type"].astype(
        pd.CategoricalDtype(categories=group_order, ordered=True)
    )
    ## Ensure correct categorical ordering by arranging var max_cell_type following group_order
    heatmap_adata = heatmap_adata[:, heatmap_adata.var_names[heatmap_adata.var["max_cell_type"].isin(group_order)]]
    heatmap_adata.var["max_cell_type"] = heatmap_adata.var["max_cell_type"].astype(
        pd.CategoricalDtype(categories=group_order, ordered=True)
    )
    heatmap_adata = heatmap_adata[:, heatmap_adata.var.sort_values("max_cell_type").index]
    ## Color
    cmap = copy.copy(cm.get_cmap("viridis"))
    cmap.set_bad(color="grey")  # NaNs will be grey
    ## Plot heatmap
    sc.pl.heatmap(
        heatmap_adata,
        var_names=heatmap_adata.var_names,
        groupby="cell_type", # Or another category for annotation
        show_gene_labels=False,
        dendrogram=False,
        figsize=(12,10),
        cmap=cmap,
        vmax=4,
        vmin=0,
        show=False
    )
    plt.savefig(f"/home/nelson.johansen/{species}_BG_heatmap.png", dpi=300)