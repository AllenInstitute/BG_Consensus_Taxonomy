import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os

import sciduck as sd

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

## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_BICAN.h5ad"))
adata.X = adata.X.astype(np.int32) ## Data is float64 from NeMO which is excessive for counts.

## Normalize the count matrix after storing the raw counts in the raw slot
adata.raw = adata.copy()

## Calculate some basic QC metrics around UMI and gene counts
sc.pp.calculate_qc_metrics(adata, inplace=True)

##
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

## Add manual QC metrics
adata.obs["total_genes"] = (adata.X > 0).sum(axis=1)

## Plot QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "total_genes"],
    jitter=0.4,
    multi_panel=True,
    save=f"_{species}_basic_qc.png",
)

## -----------------------------
## 1.1 Perform cell level basic QC
sd.qc.add_range_constraint(adata, column="total_counts", gt=2000, lt=1e6)
sd.qc.add_range_constraint(adata, column="total_genes", gt=1000, lt=15000)

## ---- FIRST SHAPE CHANGE ----
## Remove flagged cells
adata = sd.qc.apply_constraints(adata, inplace=True)
adata.obs.keeper_cells.value_counts()

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_basic_qc.h5ad"))

## -----------------------------
## 1.2 Run scVI to remove donor effects

##
adata.layers["UMIs"] = adata.raw.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")

##
adata_hvg = adata[:, adata.var.highly_variable].copy()

## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="organism_ontology_term_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
sc.tl.umap(adata, min_dist=0.3)

adata.obsm["X_umap_donor_aligned"] = adata.obsm["X_umap"].copy()

## Clean up
del adata_hvg
del adata.obsm["X_umap"]
del adata.layers["UMIs"]

## -----------------------------
## Run leiden for basic cluster-wise qc

sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=10, key_added="leiden_scVI")

## Change int64 to int32 pandas columns in adata
adata.obs = adata.obs.astype({col: 'int32' for col in adata.obs.select_dtypes(include=['int64']).columns})

## ----------------------------- 
## 1.3 cluster-wise QC metrics

## Add cluster-wise qc metrics to determine thresholds
adata.obs["cluster_doublet_score"] = adata.obs.groupby("leiden_scVI")["doublet_score"].transform("median")
adata.obs["cluster_total_counts"] = adata.obs.groupby("leiden_scVI")["total_counts"].transform("median")
adata.obs["cluster_total_genes"] = adata.obs.groupby("leiden_scVI")["total_genes"].transform("median")

## -----------------------------
## 1.3.1 Bring in old annotations from prior work
adata_AIT117 = ad.read_h5ad(os.path.join("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque", 
                                         "HMBA_Macaque_BG_082024.h5ad"), backed="r")

adata_AIT117.obs.index = adata_AIT117.obs.bc.astype(str) + "-" + adata_AIT117.obs.load_name.astype(str)
common_cells = adata_AIT117.obs.index.intersection(adata.obs.index)

adata_AIT117.obs.rename({"level4_group": "Group",
                        "level3_subclass": "Subclass",
                        "level2_class": "Class",
                        "level1_neighborhood": "Neighborhood"
    }, axis=1, inplace=True)

for anno in ["cluster", "Group", "Subclass", "Class", "Neighborhood"]:
    adata.obs["AIT117_" + anno] = "NA"
    adata.obs.loc[common_cells, "AIT117_" + anno] = adata_AIT117.obs.loc[common_cells, anno]
    adata.obs["AIT117_" + anno] = adata.obs["AIT117_" + anno].astype("category")

## Change int64 to int32 pandas columns in adata
adata.obs = adata.obs.astype({col: 'int32' for col in adata.obs.select_dtypes(include=['int64']).columns})

## 
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_cluster_qc.h5ad")) ## 
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_HMBA_cluster_qc.h5ad")) ## 

## -----------------------------
## Identify clusters that have high qc metrics but are of good quality and should be ignored when apply constraints.
## This is a manual process and should be done by eye.

exempt_clusters = ["186", "185", "170"]

## -----------------------------
##
sd.qc.add_group_level_constraint(adata, column="doublet_score", groupby = "leiden_scVI", lt=0.22, agg_func = "median")
sd.qc.add_group_level_constraint(adata, column="total_counts", groupby = "leiden_scVI", lt=2e5, agg_func = "median")
sd.qc.add_range_constraint(adata, column="doublet_score", lt=0.35)
sd.qc.add_range_constraint(adata, column="total_counts", gt=2000, lt=1e5)

adata = sd.qc.apply_constraints(adata, inplace=False)
adata.obs.keeper_cells.value_counts()

## Now we ignore the exempt clusters and remove remaining flagged cells to remove them from the dataset.

adata = adata[adata.obs["leiden_scVI"].isin(exempt_clusters) | adata.obs["keeper_cells"]].copy()

## -----------------------------
## 1.4 re-run scVI on clean-ish data
adata.layers["UMIs"] = adata.raw.X.copy() ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")

##
adata_hvg = adata[:, adata.var.highly_variable].copy()

## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="organism_ontology_term_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
sc.tl.umap(adata, min_dist=0.3)

adata.obsm["X_umap_donor_aligned"] = adata.obsm["X_umap"].copy()

## Clean up
del adata_hvg
del adata.obsm["X_umap"]
del adata.layers["UMIs"]

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_cluster_qc.h5ad")) ## Overwrites previous save.
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_HMBA_cluster_qc.h5ad")) ## 
