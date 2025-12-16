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

import os

##
species = "Macaque"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

## Save donor corrected latent space
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno_031525.h5ad"))

##
class_groups = {"CN_LGE_GABA": ["CN LGE GABA"],
                "Interneurons": ["MB Dopaminergic", "F M Glut", "NA(0.95)", "CN MGE GABA", "CN CGE GABA", "BG GABA Glut"], 
                "Non-Neurons": ["OPC-Oligo", "Astro-Epen", "Immune", "Vascular"]
                }

## Run scVI for each group
for type_label, type_list in class_groups.items():
    print(type_list)
    ### Subset to class
    adata_sub = adata[adata.obs["Class"].isin(type_list)].copy()
    adata_sub.layers["UMIs"] = adata_sub.raw.X.copy()
    ## Filter to class for both species
    sc.pp.highly_variable_genes(adata_sub, flavor="seurat_v3", n_top_genes=2000, layer="UMIs", batch_key="donor_name", subset=False)
    adata_hvg = adata_sub[:,adata_sub.var.highly_variable].copy()
    ## ----------------------------------------
    ## Integrate across donors and species
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer="UMIs",
        batch_key="donor_name",
        categorical_covariate_keys=["donor_name"],
    )
    model = scvi.model.SCVI(adata_hvg, 
                            dispersion="gene-batch", 
                            n_hidden=128, 
                            n_latent=64, 
                            n_layers=3)
    model.train(max_epochs=200)
    ## Save model
    # model.save(dir_path=wkdir + "/aibs_integration_models/")
    ## Extract integrated space
    adata_sub.obsm["X_scVI"] = model.get_latent_representation()
    ## UMAP from scVI latent space
    sc.pp.neighbors(adata_sub, use_rep="X_scVI")
    sc.tl.umap(adata_sub, min_dist=0.3)
    ## Save
    adata_sub.write_h5ad(os.path.join(cxgdir, f"{species}_basalganglia_HMBA_AIT11-9_anno_031525_{type_label}.h5ad"))
