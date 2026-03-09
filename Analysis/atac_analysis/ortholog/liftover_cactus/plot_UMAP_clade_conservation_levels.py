import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import umap
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import gzip
import glob 
import os
import anndata as ad
from scipy.sparse import csr_matrix
from ete3 import NCBITaxa
import scanpy as sc

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Load data
merged_adata = ad.read_h5ad(os.path.join(data_dir, "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"))
species_adata = merged_adata[merged_adata.obs["species"].isin(["human"])].copy()
species_adata.uns["peak_length"] = 501

## 
supplement_tables = {}

## Identify peaks which are conserved across 90% of primates (merged_adata.var.order == "Primates")
for order in species_adata.var["order"].unique():
    # if f"n_species_>=90%_aligned_in_{order}" in species_adata.obs.columns:
    #     continue
    print(f"Calculating conservation for order: {order} with {len(species_adata.var[species_adata.var['order'] == order])} species")
    species_in_order = species_adata.var[species_adata.var["order"] == order].index.tolist()
    order_adata = species_adata[:,species_in_order]
    species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] >= 0.9).sum(axis=1).A1
    species_adata.obs[f"n_species_<=10%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] <= 0.1).sum(axis=1).A1
    ## Print
    print(species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe())
    ## I need to gather the results from describe into a data.frame to save
    ## and then I can make a table for the supplement
    supplement_tables[order] = species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe()
    ## Plot
    plt.figure(figsize=(10,8))
    plt.scatter(species_adata.obsm["X_umap"][:,0], 
                species_adata.obsm["X_umap"][:,1], 
                cmap='magma',
                c=species_adata.obs[f"n_species_>=90%_aligned_in_{order}"].values,
                s=0.001, 
                alpha=0.5)
    ## Add color legend
    cbar = plt.colorbar()
    cbar.set_label(f'Number of species in {order} with >=90% alignment', rotation=270, labelpad=15)
    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.yaxis.set_ticks_position('left')
    plt.title("UMAP of CACTUS alignment to 447 species for Human peaks")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "figures", f"CACTUS_umap_{order}_>=90%_conservation.png"), dpi=900)  

for order in species_adata.var["primates_shorthand"].unique():
    # if f"n_species_>=90%_aligned_in_{order}" in species_adata.obs.columns:
    #     continue
    print(f"Calculating conservation for order: {order} with {len(species_adata.var[species_adata.var['primates_shorthand'] == order])} species")
    species_in_order = species_adata.var[species_adata.var["primates_shorthand"] == order].index.tolist()
    order_adata = species_adata[:,species_in_order]
    species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] >= 0.9).sum(axis=1).A1
    species_adata.obs[f"n_species_<=10%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] <= 0.1).sum(axis=1).A1
    ## Print
    print(species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe())
    ## I need to gather the results from describe into a data.frame to save
    ## and then I can make a table for the supplement
    supplement_tables[order] = species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe()
    ## Plot
    plt.figure(figsize=(10,8))
    plt.scatter(species_adata.obsm["X_umap"][:,0], 
                species_adata.obsm["X_umap"][:,1], 
                cmap='magma',
                c=species_adata.obs[f"n_species_>=90%_aligned_in_{order}"].values,
                s=0.001, 
                alpha=0.5)
    ## Add color legend
    cbar = plt.colorbar()
    cbar.set_label(f'Number of species in {order} with >=90% alignment', rotation=270, labelpad=15)
    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.yaxis.set_ticks_position('left')
    plt.title("UMAP of CACTUS alignment to 447 species for Human peaks")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "figures", f"CACTUS_umap_{order}_>=90%_conservation.png"), dpi=900)  

## Final conservation plots for manual primate categories
for order in species_adata.var["primates_shorthand"].cat.categories:
    # if f"n_species_>=90%_aligned_in_{order}" in species_adata.obs.columns:
    #     continue
    print(f"Calculating conservation for order: {order} with {len(species_adata.var[species_adata.var['primates_shorthand'] == order])} species")
    species_in_order = species_adata.var[species_adata.var["primates_shorthand"] == order].index.tolist()
    if f"n_species_>=90%_aligned_in_{order}" not in species_adata.obs.columns:
        order_adata = species_adata[:,species_in_order]
        species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] >= 0.9).sum(axis=1).A1
        species_adata.obs[f"n_species_<=10%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] <= 0.1).sum(axis=1).A1
    ## Print
    print(species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe())
    ## I need to gather the results from describe into a data.frame to save
    ## and then I can make a table for the supplement
    supplement_tables[order] = species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe()
    ## Plot
    plt.figure(figsize=(10,8))
    plt.scatter(species_adata.obsm["X_umap"][:,0], 
                species_adata.obsm["X_umap"][:,1], 
                cmap='Greys',
                c=species_adata.obs[f"n_species_>=90%_aligned_in_{order}"].values,
                s=0.01, 
                alpha=0.8)
    ## Add color legend
    cbar = plt.colorbar()
    cbar.set_label(f'Number of species in {order} with >=90% alignment', rotation=270, labelpad=15)
    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.yaxis.set_ticks_position('left')
    plt.title("UMAP of CACTUS alignment to 447 species for Human peaks")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "figures", f"CACTUS_umap_{order}_>=90%_conservation_final.png"), dpi=900) 

## ------------------------------------------------
## Final primates UMAP with 0 values in light blue
## ------------------------------------------------
order = "Primates"
print(f"Calculating conservation for order: {order} with {len(species_adata.var[species_adata.var['order'] == order])} species")
species_in_order = species_adata.var[species_adata.var["order"] == order].index.tolist()
order_adata = species_adata[:,species_in_order]
species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] >= 0.9).sum(axis=1).A1
species_adata.obs[f"n_species_<=10%_aligned_in_{order}"] = (order_adata.X/species_adata.uns["peak_length"] <= 0.1).sum(axis=1).A1

## Print
print(species_adata.obs[[f"n_species_>=90%_aligned_in_{order}", f"n_species_<=10%_aligned_in_{order}"]].describe())

# Define custom green to blue colormap
import matplotlib.colors as mcolors
colors = ["blue", "lightgrey", "green"]
cmap = mcolors.LinearSegmentedColormap.from_list("GreenBlue", colors)

## Plot
plt.figure(figsize=(10,8))
## Plot 0 values as light blue
# plt.scatter(species_adata.obsm["X_umap"][:,0][species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] < 40], 
#             species_adata.obsm["X_umap"][:,1][species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] < 40], 
#             color='lightblue', 
#             s=0.001, 
#             alpha=0.8)
## Plot all values as greys
plt.scatter(species_adata.obsm["X_umap"][:,0], #[species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] >= 40], 
            species_adata.obsm["X_umap"][:,1], #[species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] >= 40], 
            cmap=cmap,
            c=species_adata.obs[f"n_species_>=90%_aligned_in_{order}"].values, #[species_adata.obs[f"n_species_>=90%_aligned_in_{order}"] >= 40],
            s=0.001, 
            alpha=0.9)
## Add color legend
cbar = plt.colorbar()
cbar.set_label(f'Number of species in {order} with >=90% alignment', rotation=270, labelpad=15)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
plt.title("UMAP of CACTUS alignment to 447 species for Human peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "figures", f"CACTUS_umap_{order}_>=90%_conservation_final.png"), dpi=900)  

## Save final merged adata with conservation levels
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"
species_adata.write_h5ad(os.path.join(cxg_dir, "CACTUS", "human_zoonomia_overlap_HALPER_NCBI_all_anno_conservation_levels.h5ad"))
