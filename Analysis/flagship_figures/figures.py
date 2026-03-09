import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## Load HMBA consensus data
adata = ad.read_h5ad(os.path.join(work_dir, "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"))

n_clusters = 451
palette = sns.color_palette("hsv", n_clusters)

## Plot each species UMAP colored by cluster
for species_name in ["Macaque", "Marmoset"]:
    print(f"Processing {species_name}...")
    adata_sub = adata[adata.obs['organism'] == species_name].copy()
    ## Plot UMAP with consensus annotation colors
    colors = adata_sub.obs["Cluster"]#.map(color_map).values.astype(str)
    colors = colors.str.replace(f"{species_name}-", "", regex=True)
    colors = colors.astype(int)
    ## Create color palette for the current species
    palette = sns.color_palette("hsv", len(colors.unique()))
    ## Plot UMAP colored by cell type
    plt.figure(figsize=(12,12))
    plt.scatter(adata_sub.obsm["X_umap"][:,0], 
                adata_sub.obsm["X_umap"][:,1], 
                c=[palette[i % len(palette)] for i in colors],
                s=0.01, 
                alpha=0.7)
    plt.title(f"UMAP of HMBA BG Neurons - {species_name} - Consensus Annotation", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "flagship", "figures", f"{species_name}_hmba_neuron_umap.png"), dpi=900)

## Plot a UMAP colored Grey
plt.figure(figsize=(12,12))
plt.scatter(adata.obsm["X_umap"][:,0], 
            adata.obsm["X_umap"][:,1], 
            c="grey",
            s=0.01, 
            alpha=0.7)
plt.title(f"UMAP of HMBA BG Neurons - {species_name} - Consensus Annotation", fontsize=16)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "flagship", "figures", f"all_grey_hmba_neuron_umap.png"), dpi=900)
