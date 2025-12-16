import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import glob 
from functools import reduce
import re
import os
from collections import Counter
import scvi

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies"

##
adata = ad.read_h5ad(work_dir + "/HMBA_BG_Human_Macaque_Marmoset_alignment.h5ad", backed="r")

adata.obs.Group = adata.obs.Group.astype(str)
adata.obs.loc[adata.obs.organism == "Marmoset", "Group"] = adata.obs.loc[adata.obs.organism == "Marmoset", "Group_scanvicluster_HMBA_BG"]
adata.obs.Group = adata.obs.Group.astype('category')

##
adata_dopa = adata[adata.obs.loc[adata.obs.Group == "SN TH Dopa"].index]
adata_dopa.write_h5ad(os.path.join(work_dir, "HMBA_BG_Human_Macaque_Marmoset_alignment_dopa.h5ad"))

## Filter to class for both species
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=1000, layer="UMIs", batch_key="organism", subset=False)
adata_hvg = adata[:,adata.var.highly_variable].copy()

## ----------------------------------------
## Integrate across donors and species
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="organism", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=128, 
                        n_latent=32, 
                        n_layers=3)
model.train(max_epochs=200)
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI_xspecies"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI_xspecies")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=5, key_added="leiden_scVI")
sc.tl.umap(adata, min_dist=0.3)

##---
adata.obsm["X_umap_species_integrated"] = adata.obsm["X_umap"].copy()

##
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)


adata.write_h5ad(os.path.join(cxgdir, "HMBA_BG_Human_Macaque_Marmoset_alignment_dopa_scVI.h5ad"))

