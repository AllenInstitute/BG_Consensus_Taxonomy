import anndata as ad
import pandas as pd
import numpy as np
import umap
import os

from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/analysis/"
gini_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/specificity"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"

## Load data
merged_adata = ad.read_h5ad(os.path.join(analysis_dir, "CACTUS", "all_species_zoonomia_overlap_HALPER.h5ad"))
merged_adata.obs["region-id"] = (
    merged_adata.obs.index.astype(str) + "-" + merged_adata.obs["species"].astype(str)
)

## Load in gini scores
gini_scores = pd.read_csv(os.path.join(gini_dir, "gini_scores_combined.csv"), index_col=0)
gini_scores["region-id"] = (
    gini_scores.index.astype(str) + "-" + gini_scores["species"].astype(str)
)

##
merged_adata.obs["gini_scores"] = merged_adata.obs["region-id"].map(
    gini_scores.set_index("region-id")["gini_scores"]
)

## -----------------------------------------------------
## Plots!
## -----------------------------------------------------
plot_adata = merged_adata.copy()
print(f"Plotting UMAP with {plot_adata.n_obs} peaks")

## Plot UMAPs
species_plot = "all"

## Plot UMAP
plt.figure(figsize=(8,6))
plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c=plot_adata.obs["conservation"].map({
                "highly_conserved": "blue", 
                "actively_evolving": "orange", 
                "poorly_conserved": "red", 
                "": "lightgrey"
            }),
            s=0.0001, 
            alpha=0.5)
## Add color legend
blue_patch = mpatches.Patch(color='blue', label='highly_conserved')
orange_patch = mpatches.Patch(color='orange', label='actively_evolving')
red_patch = mpatches.Patch(color='red', label='poorly_conserved')
grey_patch = mpatches.Patch(color='lightgrey', label='other')
plt.legend(handles=[blue_patch, orange_patch, red_patch, grey_patch], markerscale=100, loc='best', fontsize='small')
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap.png"), dpi=700)

plt.figure(figsize=(8,6))
plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c = plot_adata.obs.species.map({
                "human": "blue", 
                "macaque": "orange", 
                "mouse": "red"
            }),
            s=0.0001, 
            alpha=0.5)
## Add color legend
blue_patch = mpatches.Patch(color='blue', label='human')
orange_patch = mpatches.Patch(color='orange', label='macaque')
red_patch = mpatches.Patch(color='red', label='mouse')
plt.legend(handles=[blue_patch, orange_patch, red_patch], markerscale=100, loc='best', fontsize='small')
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap_species.png"), dpi=700)

plt.figure(figsize=(8,6))
sc = plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c = plot_adata.obs["n_species_>=90%_aligned"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="#>=90%")
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap_90percent.png"), dpi=700)

plt.figure(figsize=(8,6))
sc = plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c = plot_adata.obs["n_species_<=10%_aligned"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="#<=10%")
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap_10percent.png"), dpi=700)

plt.figure(figsize=(8,6))
sc = plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c = plot_adata.obs["gini_scores"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="Gini Scores")
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap_gini_scores.png"), dpi=700)


plot_adata.obs["promoter"] = False
plot_adata.obs.loc[plot_adata.obs["human_promoter"] == 'True', "promoter"] = True
plot_adata.obs.loc[plot_adata.obs["macaque_promoter"] == 'True', "promoter"] = True
plot_adata.obs.loc[plot_adata.obs["mouse_promoter"] == 'True', "promoter"] = True


plt.figure(figsize=(8,6))
sc = plt.scatter(plot_adata.obsm["X_umap"][:, 0], 
            plot_adata.obsm["X_umap"][:, 1], 
            c = plot_adata.obs["promoter"].map({True: "blue", False: "lightgrey"}),
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="Promoter")
plt.title(f"UMAP of CACTUS alignment to 447 species for {species_plot} peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_umap_promoter.png"), dpi=700)

## -----------------------------------------------------
## KDE plots
## -----------------------------------------------------
plt.figure(figsize=(8,6))
sns.kdeplot(
    data = plot_adata.obs,
    x="n_species_<=10%_aligned", 
    y="n_species_>=90%_aligned", 
    hue="species",
    fill=False, 
    thresh=0.05, 
    levels=100,
    linewidths=0.5,
    alpha=0.7,
)
plt.title("Peak alignment KDE heatmap")
plt.xlabel("Species aligned <=10%")
plt.ylabel("Species aligned ≥90%")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"{species_plot}_peaks_kde_species_manual.png"), dpi=700)

## Look for human-specific peaks
human_specific = merged_adata[
    (merged_adata.obs["species"] == "human")
]

human_specific_strict = human_specific[
    (human_specific.obs["n_species_>=90%_aligned"] == 1) &
    (human_specific.obs["n_species_<=10%_aligned"] > 400)
].copy()