import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os
import sys
import sciduck as sd

##
species = "Human"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

os.chdir(work_dir)


adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-5_anno.h5ad")) 

df = pd.read_excel(work_dir + "/Human_AIBS_AIT19-5_consensus_anno_DateXX.xlsx", usecols=["cluster_id", "Neighborhood", "Class", 'Subclass', 'Group']) #Cluster -> cluster_id in google sheet
df.set_index("cluster_id", inplace=True)
print(df)


obs_df = adata.obs[["cluster"]].copy()
obs_df = obs_df.merge(df, left_on="cluster", right_on="cluster_id", how="left")

adata.obs["Neighborhood"] = obs_df["Neighborhood"].values
adata.obs["Class"] = obs_df["Class"].values
adata.obs["Subclass"] = obs_df["Subclass"].values
adata.obs["Group"] = obs_df["Group"].values

adata.obs["Neighborhood"] = adata.obs["Neighborhood"].str.strip() #Junk, or 'Junk '
adata.obs["Class"] = adata.obs["Class"].str.strip()
adata.obs["Subclass"] = adata.obs["Subclass"].str.strip()
adata.obs["Group"] = adata.obs["Group"].str.strip()


cols_to_replace = ['Class', 'Subclass', 'Group']
for col in cols_to_replace:
    #adata.obs[col] = adata.obs[col].cat.add_categories(['TBD']) # a few clusters' identities - 'TBD' for now
    adata.obs[col] = adata.obs[col].fillna('TBD')

sc.pl.embedding(
    adata,
    basis='X_umap',
    color=['Neighborhood', 'Class', 'Subclass', 'Group'], 
    use_raw=False,
    ncols=1
)

# Remove cells where Neighborhood is "Junk"
adata = adata[adata.obs["Neighborhood"] != "Junk"].copy() #based on Neighborhood_h5ads, more "Junk" clusters were identified.

##----------- adata.shape changed so re-run scVI.
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")
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
adata.obs = adata.obs.rename(columns={'cluster': 'cluster_id'})

adata.write_h5ad(work_dir + "/Human_basalganglia_AIBS_AIT19.5_anno_DateXXX.h5ad") 


############# Then, subset the cells from each Neighborhood_level into separate h5ad files and re-run scVI on each of them individually. 
adata_nonneuron = adata[adata.obs["Neighborhood"] == "Nonneuron"].copy()
adata_cn_lge_gaba = adata[(adata.obs["Neighborhood"] == "Subpallium GABA") & (adata.obs["Class"] == "CN LGE GABA")].copy()
adata_remaining = adata[(adata.obs["Neighborhood"] != "Nonneuron") & (adata.obs["Class"] != "CN LGE GABA")].copy()

adata_nonneuron.write_h5ad(work_dir + "/AIT19-5_Nonneuron_v2.h5ad")  # Save Group 1
adata_cn_lge_gaba.write_h5ad(work_dir + "/AIT19-5_LGE_MSNs_v2.h5ad")  # Save Group 2
adata_remaining.write_h5ad(work_dir + "/AIT19-5_Interneurons_v2.h5ad")  # Save Group 3

## Script for re-running scVI of each Neighborhood_level h5ad was omitted here, as it is identical to the latter part of Step 7. For more details, please see Step 7.


