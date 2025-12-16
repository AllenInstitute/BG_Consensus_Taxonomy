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

df = pd.read_excel(work_dir + "/Human_AIBS_AIT19-5_consensus_anno_table.xlsx", usecols=["Cluster", "Neighborhood", "Class"]) #This .xlsx file was the google sheet. Using this, we did manual annotation per species. Now, all neighborhoold level tags were assigned based on mapping results and markers. file was also downloaded into my home dir.
df.set_index("Cluster", inplace=True)
print(df)

obs_df = adata.obs[["cluster"]].copy()
obs_df = obs_df.merge(df, left_on="cluster", right_on="Cluster", how="left")

adata.obs["Neighborhood"] = obs_df["Neighborhood"].values
adata.obs["Class"] = obs_df["Class"].values
adata.obs["Neighborhood"] = adata.obs["Neighborhood"].str.strip() #Junk, or 'Junk '
adata.obs["Class"] = adata.obs["Class"].str.strip()



# Remove cells where Neighborhood is "Junk"
adata = adata[adata.obs["Neighborhood"] != "Junk"].copy() #'Junk' were assigned when going through the annotation table (QC metrics + mapping results + marker expresssion + species alignment)

# Define three groups based on Neighborhood and Class
# Group 1: Cells where Neighborhood == "Nonneuron"
adata_nonneuron = adata[adata.obs["Neighborhood"] == "Nonneuron"].copy()

# Group 2: Cells where Neighborhood == "Subpallium GABA" AND Class == "CN LGE GABA"
adata_cn_lge_gaba = adata[
    (adata.obs["Neighborhood"] == "Subpallium GABA") & 
    (adata.obs["Class"] == "CN LGE GABA")
].copy()

# Group 3: All remaining cells (i.e., not in Group 1 and not in Group 2)
adata_remaining = adata[
    (adata.obs["Neighborhood"] != "Nonneuron") & 
    (adata.obs["Class"] != "CN LGE GABA")
].copy()

# Save each group as a separate h5ad file** (human BG with over 1M cells - too big to run in a same slurm job.
adata_nonneuron.write_h5ad(work_dir + "/AIT19-5_Nonneuron.h5ad")  # Save Group 1
adata_cn_lge_gaba.write_h5ad(work_dir + "/AIT19-5_LGE_MSNs.h5ad")  # Save Group 2
adata_remaining.write_h5ad(work_dir + "/AIT19-5_Interneurons.h5ad")  # Save Group 3

######---------------------------------- re-run scVI ----------------------------------######
###------- (1) re-run scVI for Interneurons of AIT19-5. Then do carefully manual annotations based on this new h5ad.

#adata = sc.read_h5ad(work_dir + "/AIT19-5_Interneurons.h5ad") 
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=2000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")
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

adata.write_h5ad(work_dir + "/AIT19-5_Interneurons_scVI.h5ad") 

###------- (2) re-run scVI for MSNs+Granular cells of AIT19-5.

#adata = sc.read_h5ad(work_dir + "/AIT19-5_LGE_MSNs.h5ad") 
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
adata = adata[adata.obs["donor_id"] != "H19.30.004"].copy() #only 2 cells in MSNs for this donor (not sampled from STR).
print(adata.shape)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")
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

adata.write_h5ad(work_dir + "/AIT19-5_LGE_MSNs_scVI.h5ad") 

###------- (3) re-run scVI for NNs
#adata = sc.read_h5ad(work_dir + "/AIT19-5_Nonneuron.h5ad") 
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=2000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")
##
adata_hvg = adata[:, adata.var.highly_variable].copy()

scvi.model.SCVI.setup_anndata(adata_hvg, layer="UMIs", batch_key="donor_id", categorical_covariate_keys=["donor_id"])

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene", #cannot set the "gene-batch"
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

adata.write_h5ad(work_dir + "/AIT19-5_Nonneuron_scVI.h5ad") 

#### then, Use these new h5ad files — referencing their clustering patterns on UMAP and the expression levels of marker genes —along with the mapping results from the annotation file to manually annotate each cluster one by one.
#### Keep updating the "Human_AIBS_AIT19-5_consensus_anno_table.xlsx" annotation table.
