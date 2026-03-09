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
import scvi

from tqdm import tqdm

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"


## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
# adata = ad.read_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_neuron_alignment_labeltransfer_allgenes.h5ad"), backed="r")

## -------------------------
# group_filter = ["STR-BF TAC3-PLPP4-LHX8 GABA", "STR TAC3-PLPP4 GABA"]
# mouse_filter = ["0840 STR Lhx8 Gaba_1", "0841 STR Lhx8 Gaba_1", "0842 STR Lhx8 Gaba_1"]

# adata_tac3 = adata[adata.obs["Group"].isin(group_filter) | adata.obs["AIT21.cluster"].isin(mouse_filter)]
# adata_tac3.write_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_TAC3_alignment_labeltransfer_allgenes.h5ad"))

##
adata = ad.read_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_TAC3_alignment_labeltransfer_allgenes.h5ad"))

## Normalize and log
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

## neighbor graph
# sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
# sc.tl.umap(adata, min_dist=0.2)

##
adata.write_h5ad(os.path.join(cxg_dir, "BasalGanglia", "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_TAC3_alignment_labeltransfer_allgenes.h5ad"))

##
metadata = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_TAC3_alignment_labeltransfer_allgenes_annotations/anno.csv", skiprows=2, index_col=0)
metadata = metadata.loc[metadata.keepers == True, :]

adata = adata[metadata.index, :].copy()

## Fix labels
adata.obs.loc[adata.obs["Cluster"] == "Macaque-124", "Group"] = "STR-BF TAC3-PLPP4-LHX8 GABA"

## -------------------------------
## Plot UMAP colored by cell type
plt.figure(figsize=(12,12))
adata_plot = adata[adata.obs["organism"] != "Mouse", :]
plt.scatter(adata_plot.obsm["X_umap"][:,0], 
            adata_plot.obsm["X_umap"][:,1], 
            c=adata_plot.obs["Group"].map(color_map).values.astype(str),
            s=5, 
            alpha=0.8)
adata_plot = adata[adata.obs["organism"] == "Mouse", :]
plt.scatter(adata_plot.obsm["X_umap"][:,0], 
            adata_plot.obsm["X_umap"][:,1], 
            c="#fc751c",
            s=5, 
            alpha=0.8)
## Add color legend
plt.title("UMAP of HMBA BG Neurons - Consensus Annotation", fontsize=16)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "vignette", "TAC3", f"xpsecies_TAC3_umap.png"), dpi=900)  
plt.close()

## ---------------------------------
## Homology
adata_homology = adata.copy()
adata_homology.obs["Group"] = adata_homology.obs["Group"].astype(str)
adata_homology.obs.loc[adata_homology.obs["Cluster"].isin(["Macaque-119", "Macaque-118", "Macaque-122", "Marmoset-283", "Human-183"]), "Group"] = "STR-BF TAC3-PLPP4-LHX8 GABA_2"

## Define the annotation level from the "reference"
category_set_ref = adata_homology.obs["Group"].dropna().unique().tolist()
category_set_ref = [x for x in category_set_ref if x != 'nan']

## Define the annotation level from the "query"
category_set_query = adata_homology.obs["AIT21.cluster"].dropna().unique().tolist()
category_set_query = [x for x in category_set_query if x != 'nan']

## Another way to define category set when you want to include only clusters with > 100 cells for example.
# counts = adata_homology.obs.cluster_id.value_counts()
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
for anno in adata_homology.obs["Group"].dropna().unique().tolist():
    homology_df["ref"]["centroid"][anno] = np.median(adata_homology[(adata_homology.obs["Group"] == anno) & (adata_homology.obs.organism != "Mouse")].obsm["X_scVI"], axis=0)  

## Compute centroids for query.
for anno in adata_homology.obs["AIT21.cluster"].dropna().unique():
    homology_df["query"]["centroid"][anno] = np.median(adata_homology[(adata_homology.obs["AIT21.cluster"] == anno) & (adata_homology.obs.organism == "Mouse")].obsm["X_scVI"], axis=0)

## Build homology matrix
for anno in tqdm(category_set_ref):
    for anno_query in adata_homology.obs["AIT21.cluster"].dropna().unique():
        if np.isnan(homology_df["ref"]["centroid"][anno]).any() | np.isnan(homology_df["query"]["centroid"][anno_query]).any():
            homology_df["homology"][anno][anno_query] = 0
        else:
            homology_df["homology"][anno][anno_query] = scipy.spatial.distance.euclidean(homology_df["ref"]["centroid"][anno], homology_df["query"]["centroid"][anno_query])

##
homology_mat = homology_df["homology"].loc[category_set_query, category_set_ref]

## min-max normalization
homology_mat_norm = 1 - (homology_mat - homology_mat.min().min()) / ((homology_mat.max().max() - homology_mat.min().min()) + 1e-9)
homology_mat_norm = homology_mat_norm.fillna(0)

##
row_best_match = homology_mat_norm.idxmax(axis=1)

# Sort rows based on the index of that best-matching column in your predefined order
column_order = cluster_meta.Group[cluster_meta.Group.isin(homology_mat_norm.columns).tolist()].tolist() + ["STR-BF TAC3-PLPP4-LHX8 GABA_2"]
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
plt.savefig(os.path.join(analysis_dir, "vignette", "TAC3", "homology_hmba_bg_mouse_TAC3_diagonalized_fixedcols.png"), dpi=900)
plt.savefig(os.path.join(analysis_dir, "vignette", "TAC3", "homology_hmba_bg_mouse_TAC3_diagonalized_fixedcols.pdf"), dpi=900)
plt.close()

## -------------------------------
## Cluster level expression of marker genes
adata.obs["cluster_plot"] = adata.obs["Cluster"].astype(str)
adata.obs.loc[adata.obs["cluster_plot"] == "nan", "cluster_plot"] = adata.obs.loc[adata.obs["cluster_plot"] == "nan", "AIT21.cluster"]

marker_genes = {
    "broad": ["LHX8", "TAC3", "PLPP4"],
    "_2": ["SOX5", "SORCS1"],
    "_1": ["SOX6", "IGF1", "GLP1R", "RUNX1T1"]
}

markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)