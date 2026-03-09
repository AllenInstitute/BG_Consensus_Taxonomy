import anndata as ad
import pandas as pd
import numpy as np
import umap
import os

from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

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

## -----------------------------------------------------
## Species-specific UMAP
## -----------------------------------------------------
species_color_map = {
    "human": "blue",
    "macaque": "green",
    "mouse": "orange",
}

for species in merged_adata.obs["species"].unique():
    ## colors
    colors = np.where(merged_adata.obs["species"].values == species, species_color_map[species], "lightgrey")
    ## Plot
    plt.figure(figsize=(6, 6))
    plt.scatter(
        merged_adata.obsm["X_umap"][:, 0],
        merged_adata.obsm["X_umap"][:, 1],
        c=colors,
        s=0.0001,  
        alpha=0.5,
        linewidths=0, # no edge
        rasterized=True
    )
    ## Legend (proxy patches)
    species_patch = mpatches.Patch(color=species_color_map[species], label=species)
    grey_patch = mpatches.Patch(color='lightgrey', label='other species')
    plt.legend(handles=[species_patch, grey_patch], loc='best', fontsize='small')
    ## Labels and save
    plt.title("UMAP of CACTUS alignment to 447 species for all peaks")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    ##
    plt.savefig(os.path.join(analysis_dir, "figures", f"CACTUS_umap_{species}.png"), dpi=700)  

## -----------------------------------------------------
## Actively evolving >=30% <=70%
## -----------------------------------------------------
merged_adata.uns["peak_length"] = 501
merged_adata.obs["n_species_>=60%_aligned"] = (merged_adata.X/merged_adata.uns["peak_length"] >= 0.6).sum(axis=1).A1
merged_adata.obs["n_species_<=80%_aligned"] = (merged_adata.X/merged_adata.uns["peak_length"] <= 0.8).sum(axis=1).A1

merged_adata.obs["actively_evolving"] = (
    (merged_adata.obs["n_species_>=60%_aligned"] >= 1) &
    (merged_adata.obs["n_species_<=80%_aligned"] >= 1)
)

plt.figure(figsize=(8,6))

mask = merged_adata.obs["actively_evolving"] == False
sc = plt.scatter(merged_adata.obsm["X_umap"][mask, 0], 
            merged_adata.obsm["X_umap"][mask, 1], 
            color="lightgrey",
            s=0.0001,
            alpha=0.5)

mask = merged_adata.obs["actively_evolving"] == True
sc = plt.scatter(merged_adata.obsm["X_umap"][mask, 0], 
            merged_adata.obsm["X_umap"][mask, 1], 
            color="red",
            s=0.0001,
            alpha=0.5)

## Add color legend
plt.colorbar(sc, label="#>=30%")
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_actively_evolving.png"), dpi=700)

## -----------------------------------------------------
## >=90%, <=10%, evo distance
## -----------------------------------------------------
plt.figure(figsize=(8,6))
sc = plt.scatter(merged_adata.obsm["X_umap"][:, 0], 
            merged_adata.obsm["X_umap"][:, 1], 
            c = merged_adata.obs["n_species_>=90%_aligned"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="#>=90%")
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_90percent.png"), dpi=700)

plt.figure(figsize=(8,6))
sc = plt.scatter(merged_adata.obsm["X_umap"][:, 0], 
            merged_adata.obsm["X_umap"][:, 1], 
            c = merged_adata.obs["n_species_<=10%_aligned"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="#<=10%")
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_10percent.png"), dpi=700)

plt.figure(figsize=(8,6))
sc = plt.scatter(merged_adata.obsm["X_umap"][:, 0], 
            merged_adata.obsm["X_umap"][:, 1], 
            c = merged_adata.obs["evo_distance"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="#>=90%")
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_evo_distance.png"), dpi=700)

## -----------------------------------------------------
## Gini scores and cell type calls
## -----------------------------------------------------
plt.figure(figsize=(8,8))

## Plot non-TE first (grey)
non_marker_mask = merged_adata.obs["gini_scores"] <= 0.8
plt.scatter(
    merged_adata.obsm["X_umap"][non_marker_mask, 0],
    merged_adata.obsm["X_umap"][non_marker_mask, 1],
    color="lightblue",
    s=0.1,
    alpha=0.5,
    label="Non-Markers",
    edgecolors='none'
)

non_marker_mask = merged_adata.obs["promoter"] == 'True'
non_marker_mask = non_marker_mask & (merged_adata.obs["species"] == "human")
plt.scatter(
    merged_adata.obsm["X_umap"][non_marker_mask, 0],
    merged_adata.obsm["X_umap"][non_marker_mask, 1],
    color="blue",
    s=0.25,
    alpha=0.6,
    label="Promoters",
    edgecolors='none'
)

## Plot TE on top (red)
marker_mask = merged_adata.obs["gini_scores"] > 0.8
marker_mask = marker_mask & (merged_adata.obs["species"] == "human")
print(marker_mask.sum())
plt.scatter(
    merged_adata.obsm["X_umap"][marker_mask, 0],
    merged_adata.obsm["X_umap"][marker_mask, 1],
    color="red",
    s=1.0,
    alpha=0.6,
    label="Markers",
    edgecolors='none'
)

## Plot TE on top (red)
marker_mask = merged_adata.obs["gini_scores"] > 0.9
marker_mask = marker_mask & (merged_adata.obs["species"] == "human")
print(marker_mask.sum())
plt.scatter(
    merged_adata.obsm["X_umap"][marker_mask, 0],
    merged_adata.obsm["X_umap"][marker_mask, 1],
    color="darkred",
    s=1.0,
    alpha=0.6,
    label="Markers",
    marker='o',
    edgecolors='none'
)

## Add color legend
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_gini_scores.png"), dpi=900)

## -----------------------------------------------------
## Plot top 100 peaks per cell type per species
## -----------------------------------------------------

## Load in region and cell type metadata
top_peaks = pd.read_csv(os.path.join("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/specificity", 
                                     "top_100_peaks_per_cell_type_per_species.csv"))

## Merge into adata obs
plot_adata = merged_adata[merged_adata.obs_names.isin(top_peaks["region"])]
plot_adata.obs["is_top_peak"] = 'True'
plot_adata.obs["Group"] = plot_adata.obs_names.map(
    top_peaks.set_index("region")["cell_type"]
)

## Annotation up the taxonomy from Group to Subclass and Class
plot_adata.obs["Subclass"] = plot_adata.obs["Group"].map(
    cluster_meta.set_index("Group")["Subclass"]
)
plot_adata.obs["Class"] = plot_adata.obs["Group"].map(
    cluster_meta.set_index("Group")["Class"]
)

## Plot TE on top (red)
plt.figure(figsize=(8,8))

# Base scatter (all peaks, colored by Class)
class_colors = {
    "GABAergic": "blue",
    "Glutamatergic": "orange",
    "Non-Neurons": "green",
    "Motor Neurons": "red",
    "Cholinergic": "purple"
}
colors = plot_adata.obs["Class"].map(class_colors).values
colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values

plt.scatter(
    plot_adata.obsm["X_umap"][:, 0],
    plot_adata.obsm["X_umap"][:, 1],
    c=colors,
    s=3,
    alpha=0.6,
    marker='o',
    edgecolors='none',
    rasterized=True
)

## Add color legend
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_top_peaks.png"), dpi=900)



## -----------------------------------------------------
## Not useful
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

plt.figure(figsize=(8, 8))

# Base color mapping
class_colors = {
    "GABAergic": "blue",
    "Glutamatergic": "orange",
    "Non-Neurons": "green",
    "Motor Neurons": "red",
    "Cholinergic": "purple"
}

# Create DataFrame with UMAP coordinates and class labels
umap_df = pd.DataFrame({
    "UMAP1": plot_adata.obsm["X_umap"][:, 0],
    "UMAP2": plot_adata.obsm["X_umap"][:, 1],
    "Class": plot_adata.obs["Class"].values
})

# Plot KDE for each class
for cls, color in class_colors.items():
    plt.figure(figsize=(8, 8))
    cls_df = umap_df[umap_df["Class"] == cls]
    if len(cls_df) > 10:  # skip very small groups
        sns.kdeplot(
            data=cls_df,
            x="UMAP1",
            y="UMAP2",
            fill=False,
            thresh=0.05,
            color=color,
            linewidths=1.2,
            alpha=0.8,
            levels=100
        )
    # Legend
    handles = [
        plt.Line2D([], [], color=c, marker='s', linestyle='None', markersize=8, label=cls)
        for cls, c in class_colors.items()
    ]
    plt.legend(handles=handles, loc='best', fontsize='small', frameon=False)
    # Labels and Save
    plt.title("UMAP of CACTUS alignment to 447 species for all peaks (KDE density)")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_top_peaks_KDE_{cls}.png"), dpi=900)
    plt.show()


## -----------------------------------------------------
## Promoter, phyloP, TEs
## -----------------------------------------------------
plt.figure(figsize=(8,8))
sc = plt.scatter(merged_adata.obsm["X_umap"][:, 0], 
            merged_adata.obsm["X_umap"][:, 1], 
            c = merged_adata.obs["promoter"].map({'True': "blue", 'False': "lightgrey", 'nan': "lightgrey"}),
            s=0.0001, 
            alpha=0.5)
## Add color legend
blue_patch = mpatches.Patch(color='blue', label='promoter')
grey_patch = mpatches.Patch(color='lightgrey', label='non-promoter')
plt.legend(handles=[blue_patch, grey_patch], markerscale=100, loc='best', fontsize='small')
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_promoter.png"), dpi=700)

##
plt.figure(figsize=(8,6))
sc = plt.scatter(merged_adata.obsm["X_umap"][:, 0], 
            merged_adata.obsm["X_umap"][:, 1], 
            c = merged_adata.obs["phyloP_mean"],
            cmap = "magma",
            s=0.0001, 
            alpha=0.5)
## Add color legend
plt.colorbar(sc, label="phyloP_mean")
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_phyloP_mean.png"), dpi=700)

## ------------------------------------------------------
## Plot TE
## ----------------------------------------------------
plt.figure(figsize=(6,6))

## Plot non-TE first (grey)
non_te_mask = merged_adata.obs["TE"] != 'True'
plt.scatter(
    merged_adata.obsm["X_umap"][non_te_mask, 0],
    merged_adata.obsm["X_umap"][non_te_mask, 1],
    color="blue",
    s=0.0001,
    alpha=0.3,
    label="non-TE"
)

## Plot TE on top (red)
te_mask = merged_adata.obs["TE"] == 'True'
plt.scatter(
    merged_adata.obsm["X_umap"][te_mask, 0],
    merged_adata.obsm["X_umap"][te_mask, 1],
    color="red",
    s=0.0001,
    alpha=0.6,
    label="TE"
)

## Legend
red_patch = mpatches.Patch(color='red', label='TE')
grey_patch = mpatches.Patch(color='lightgrey', label='non-TE')
plt.legend(handles=[red_patch, grey_patch], markerscale=100, loc='best', fontsize='small')

## Labels and save
plt.title(f"UMAP of CACTUS alignment to 447 species for all peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_umap_TE.png"), dpi=700)












## -----------------------------------------------------
## KDE plots
## -----------------------------------------------------
plt.figure(figsize=(8,6))
sns.kdeplot(
    data = merged_adata.obs,
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
plt.savefig(os.path.join(analysis_dir, "CACTUS", "figures", f"CACTUS_peaks_kde_species_manual.png"), dpi=700)

## Look for human-specific peaks
human_specific = merged_adata[
    (merged_adata.obs["species"] == "human")
]

human_specific_strict = human_specific[
    (human_specific.obs["n_species_>=90%_aligned"] == 1) &
    (human_specific.obs["n_species_<=10%_aligned"] > 400)
].copy()