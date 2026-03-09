import os
import anndata as ad
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
import pertpy as pt
import scanpy as sc
import scipy

from tqdm import tqdm

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
adata = ad.read_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_neuron_alignment_labeltransfer_allgenes.h5ad"), backed="r")

## -------------------------

## Define the annotation level from the "reference"
category_set_ref = adata.obs["Group"].dropna().unique().tolist()
category_set_ref = [x for x in category_set_ref if x != 'nan']

## Define the annotation level from the "query"
category_set_query = adata.obs["AIT21.supertype"].dropna().unique().tolist()
category_set_query = [x for x in category_set_query if x != 'nan']

high_count_query = adata.obs["AIT21.supertype"].value_counts()[lambda x: x > 100].index.tolist()

## Another way to define category set when you want to include only clusters with > 100 cells for example.
# counts = adata.obs.cluster_id.value_counts()
# category_set_query = counts[counts >= 100].index.tolist()

## Full set of categories
category_set = np.unique(category_set_ref + category_set_query)

## Get centroids for each anno in scVI integrated space
homology_df = {
    "ref": pd.DataFrame(index=category_set_ref, columns=["centroid"]),
    "query": pd.DataFrame(index=category_set_query, columns=["centroid"]),
    "homology": pd.DataFrame(index=category_set, 
                                columns=category_set)
}

## Compute centroids for reference.
for anno in adata.obs["Group"].dropna().unique().tolist():
    homology_df["ref"]["centroid"][anno] = np.median(adata[(adata.obs["Group"] == anno) & (adata.obs.organism != "Mouse")].obsm["X_scVI"], axis=0)  

## Compute centroids for query.
for anno in adata.obs["AIT21.supertype"].dropna().unique():
    homology_df["query"]["centroid"][anno] = np.median(adata[(adata.obs["AIT21.supertype"] == anno) & (adata.obs.organism == "Mouse")].obsm["X_scVI"], axis=0)

## Build homology matrix
for anno in tqdm(category_set_ref):
    for anno_query in adata.obs["AIT21.supertype"].dropna().unique():
        if np.isnan(homology_df["ref"]["centroid"][anno]).any() | np.isnan(homology_df["query"]["centroid"][anno_query]).any():
            homology_df["homology"][anno][anno_query] = 0
        else:
            homology_df["homology"][anno][anno_query] = scipy.spatial.distance.euclidean(homology_df["ref"]["centroid"][anno], homology_df["query"]["centroid"][anno_query])

##
homology_mat = homology_df["homology"].loc[category_set_query, category_set_ref]

## min-max normalization
homology_mat_norm = 1 - (homology_mat - homology_mat.min().min()) / ((homology_mat.max().max() - homology_mat.min().min()) + 1e-9)
homology_mat_norm = homology_mat_norm.fillna(0)
homology_mat_norm = homology_mat_norm.loc[high_count_query, :]

##
row_best_match = homology_mat_norm.idxmax(axis=1)

# Sort rows based on the index of that best-matching column in your predefined order
column_order = cluster_meta.Group[cluster_meta.Group.isin(homology_mat_norm.columns).tolist()]
row_order = row_best_match.map({col: i for i, col in enumerate(column_order)}).sort_values().index

## Apply both orderings
homology_mat_ordered = homology_mat_norm.loc[row_order, column_order]

## 
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
colors = LinearSegmentedColormap.from_list("white_black", ["white", "black"])

## plot heatmap
plt.figure(figsize=(20, 20))
sns.heatmap(
    homology_mat_ordered,
    cmap=colors,
    norm=PowerNorm(gamma=3.5),
    vmin=0,
    vmax=1,
    xticklabels=True,
    yticklabels=True,
    #square=True,
    linewidths=0.3,
    linecolor="white",
    cbar_kws={"label": "1 - normalized homology"},
)
plt.title("Homology Matrix (Aligned to Column Order)", fontsize=16)
plt.xlabel("HMBA BG Neuron Groups", fontsize=12)
plt.ylabel("Mouse BG Neuron Supertypes (aligned)", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "homology", "figures", "homology_hmba_bg_mouse_diagonalized_fixedcols.png"), dpi=900)
plt.savefig(os.path.join(analysis_dir, "homology", "figures", "homology_hmba_bg_mouse_diagonalized_fixedcols.pdf"), dpi=900)
plt.close()

## Save to plot in R.
homology_df["homology"].to_csv("homology_hmba_bg_mouse.csv")