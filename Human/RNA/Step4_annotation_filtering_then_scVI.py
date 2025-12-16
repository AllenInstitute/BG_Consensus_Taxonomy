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

species = "Human"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

os.chdir(work_dir)


adata = sc.read("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies/HMBA_BG_Human_Macaque_alignment.h5ad")
adata

#leiden_align_to_remove_list = ["2", "3", "4", "13", "14", "16", "17", "19", "20", "22", "23", "24", "36", "42", "43", "57", "91", "94", "95", "108", "113", "115", "119", "120", "126", "127", "128",
#                               "129", "130", "131", "133", "134", "136", "140", "141", "143"] #https://docs.google.com/spreadsheets/d/1TNgwW_tw67Kq6kt-mxRFPVyo38Oqba00ZtlZpSA2TRw/edit?gid=0#gid=0

leiden_align_to_remove_list = ["2", "3", "4", "13", "14", "16", "17", "20", "22", "23", "24", "36", "42", "43", "57", "91", "108", "119", "120", "126", "127", "129", "130", "131", "133", "136", "140", "143"]

adata.obs["leiden_align_status"] = adata.obs["leiden_scVI"].apply(lambda x: "remove" if x in leiden_align_to_remove_list else "keep")
adata.obs["leiden_align_status"].value_counts() #same with awk '{sum += $N} END {print sum}' for google sheet # 180264

df = adata.obs[adata.obs['species'] == "Homo Sapiens"]
print(df.shape)
df.to_csv(os.path.join(work_dir, f"leiden_align_status_for_Human.csv"))

del adata

adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno.h5ad"))
print(adata.shape)

assert set(df.index) == set(adata.obs.index), "df index does not match adata.obs"
adata.obs.loc[df.index, "leiden_align_status"] = df["leiden_align_status"]
print(adata.obs["leiden_align_status"].head())


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


sc.pl.embedding(adata, basis='X_umap’, color=['leiden_align_remove_ratio_per_cluster’],  use_raw=False, ncols=1)

plt.figure(figsize=(8, 5))
# Histogram of remove_ratio values
sns.histplot(remove_ratios["remove_ratio"], bins=20, kde=True, color="steelblue")
plt.xlabel("Remove Ratio")
plt.ylabel("Number of Clusters")
plt.title("Distribution of Remove Ratios Across Clusters")
plt.show()


adata = adata[
    (adata.obs["leiden_align_remove_ratio_per_cluster"] <= 0.73) & 
    ((adata.obs["AIT193_MapMyCells_Neighborhood_label"] == "Nonneuron") | 
    ((adata.obs["AIT193_MapMyCells_Neighborhood_label"] != "Nonneuron") & (adata.obs["gene_counts_0"] >= 2000)))
].copy()
print(adata.shape)


sc.pl.embedding(
    adata,
    basis='X_umap',
    color=['leiden_align_remove_ratio_per_cluster', 'AIT193_MapMyCells_Neighborhood_label'], 
    use_raw=False,
    ncols=2
)

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno_filtered.h5ad"))



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

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno_filtered_scVI.h5ad"))

