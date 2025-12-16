import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os
import sys
import sciduck as sd
import matplotlib.pyplot as plt
import seaborn as sns


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

##
adata = sc.read("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies/HMBA_BG_Human_Macaque_alignment.h5ad", backed="r")
adata

align_qc_sheet = pd.read_csv(os.path.join(work_dir, "Human_Macaque_alignment_qc.csv"))
leiden_align_to_remove_list = align_qc_sheet.loc[align_qc_sheet["to_remove"] == "YES", "leiden_scVI"].astype(str).tolist()

adata.obs["leiden_align_status"] = adata.obs["leiden_scVI"].apply(lambda x: "remove" if x in leiden_align_to_remove_list else "keep")
adata.obs["leiden_align_status"].value_counts() #same with awk '{sum += $N} END {print sum}' for google sheet # 180264

df = adata.obs[adata.obs['organism'] == species]
print(df.shape)
df.to_csv(os.path.join(work_dir, f"leiden_align_status_for_Macaque.csv"))

del adata

## ------------------------------ 
## 
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-8_anno.h5ad"))
print(adata.shape)

alignment_results = pd.read_csv(os.path.join(work_dir, f"leiden_align_status_for_Macaque.csv"), index_col=0)
assert set(alignment_results.index) == set(adata.obs.index), "df index does not match adata.obs"
adata.obs.loc[alignment_results.index, "leiden_align_status"] = alignment_results["leiden_align_status"]
adata.obs.loc[alignment_results.index, "leiden_align_scVI"] = alignment_results["leiden_scVI"]

##  
cluster_counts = adata.obs.groupby("cluster").size() ## total number of cells each cluster
remove_counts = adata.obs.groupby("cluster")["leiden_align_status"].apply(lambda x: (x == "remove").sum()) ## the number of "remove" cells each cluster

##  Calculate the proportion of "remove" cells in each cluster
remove_ratios = (remove_counts / cluster_counts).reset_index()
remove_ratios.columns = ["cluster", "remove_ratio"]
remove_ratios.to_csv(os.path.join(work_dir, f"leiden_align_remove_ratios_per_cluster.csv"))

##  Map the remove_ratio values back to adata.obs based on the 'cluster'
adata.obs["leiden_align_remove_ratio_per_cluster"] = adata.obs["cluster"].map(remove_ratios.set_index("cluster")["remove_ratio"])

##  Verify that the new column is correctly added
print(adata.obs[["cluster", "leiden_align_remove_ratio_per_cluster"]].head())

##
sc.pl.embedding(adata, basis="X_umap_donor_aligned", color=["leiden_align_remove_ratio_per_cluster"],  use_raw=False, ncols=1, save=f"{species}_leiden_align_remove_ratio_per_cluster.png")

plt.figure(figsize=(8, 5))
# Histogram of remove_ratio values
sns.histplot(remove_ratios["remove_ratio"], bins=20, kde=True, color="steelblue")
plt.xlabel("Remove Ratio")
plt.ylabel("Number of Clusters")
plt.title("Distribution of Remove Ratios Across Clusters")
plt.savefig(os.path.join("figures", f"{species}_leiden_align_remove_ratio_per_cluster_hist.png"))

##
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT11-8_anno_qc.h5ad"))


adata.obs.leiden_align_scVI = adata.obs.leiden_align_scVI.astype(int).astype(str).astype('category')

# Assume adata.obs contains your metadata
cluster_col = "AIT117_MapMyCells_Group_label"  # Replace with actual cluster column name
bool_col = "leiden_align_scVI"  # Replace with actual True/False column name

# Convert to DataFrame
df = adata.obs[[bool_col, cluster_col]]
df = df.loc[df.leiden_align_scVI.isin(leiden_align_to_remove_list)]
df.leiden_align_scVI = df.leiden_align_scVI.cat.remove_unused_categories()

# Count occurrences of True/False per cluster
counts = df.groupby([bool_col, cluster_col]).size().unstack(fill_value=0)

# Normalize to get proportions
proportions = counts.div(counts.sum(axis=1), axis=0)

# Create heatmap
plt.figure(figsize=(20, 10))
sns.heatmap(proportions, annot=False, cmap="coolwarm", fmt=".2f", linewidths=0.5)

# Customizing the plot
plt.xlabel("Cluster")
plt.ylabel("Proportion")
plt.title("Proportion of True/False in Each Cluster")
plt.legend(title=bool_col, bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(rotation=90, ha="right") 
plt.tight_layout()
plt.savefig(os.path.join("figures", f"{species}_leiden_scVI_per_group.png"))

##
adata = adata[
    (adata.obs["leiden_align_remove_ratio_per_cluster"] <= 0.7) & 
    ((adata.obs["AIT117_MapMyCells_Neighborhood_label"] == "Nonneuron") | 
    ((adata.obs["AIT117_MapMyCells_Neighborhood_label"] != "Nonneuron") & (adata.obs["total_genes"] >= 2000)))
].copy()
print(adata.shape)


sc.pl.embedding(
    adata,
    basis='X_umap',
    color=['leiden_align_remove_ratio_per_cluster', 'AIT117_MapMyCells_Neighborhood_label'], 
    use_raw=False,
    ncols=2
)

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT11-8_anno_filtered.h5ad"))



######---------------------------------- re-run scVI ----------------------------------######

adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")

##
adata_hvg = adata[:, adata.var.highly_variable].copy()

## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="donor_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=128, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()

del adata_hvg
del adata.layers["UMIs"]


## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata = sc.tl.umap(adata, min_dist=0.3, copy=True)

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT11-8_anno_filtered_scVI.h5ad"))
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT11-8_anno_filtered_scVI.h5ad"))



# ## Load expression data
# adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-8_anno.h5ad"))

# ## Lets gather columns to put into annotation table
# anno = adata.obs.columns[adata.obs.columns.str.contains("|".join(["AIT117_.*_label", "AIT193_.*_label", "Siletti_.*_label", "ABCmouse_.*_label"]), regex=True)].tolist()
# numeric_anno = adata.obs.columns[adata.obs.columns.str.contains("|".join([".*_avg_correlation", ".*_bootstrapping_probability"]), regex=True)].tolist()

# ## Previous work
# anno_prev = ["AIT117_Neighborhood", "AIT117_Class", "AIT117_Subclass", "AIT117_Group", "AIT117_cluster"]

# ## Additional metadata
# meta = ["donor_name", "anatomical_region", "sex", "age", "barcoded_cell_sample_label"]

# ##
# anno_table = sd.anno.build_annotation_table(adata, 
#                                     group_by="cluster", 
#                                     categorical_annotations=anno_prev + anno + meta,
#                                     numeric_annotations=["doublet_score", "total_genes", "total_counts", "percent_keeper"] + numeric_anno, 
#                                     min_percent=0.05, 
#                                     annotation_alerts={"donor_name": 0.90})
# anno_table.to_csv(os.path.join(work_dir, f"{species}_{pipeline}_consensus_anno_table.csv"))

##
