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
import matplotlib.pyplot as plt

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))


## Load HMBA consensus data
adata = ad.read_h5ad(os.path.join(work_dir, "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"))

## 
grouping = {
    "LGE": ["CN LGE GABA"],
    "IN": ["F M GABA", "Cx GABA", "CN GABA-Glut", "CN CGE GABA", "CN MGE GABA", "F M Glut", "M Dopa"],
    "NN": ["OPC-Oligo", "Astro-Epen", "Immune", "Vascular"]
}

## Subset to each Class and run HVG, scVI, UMAP
for group in grouping.keys():
    print(f"Processing {group}...")
    adata_sub = adata[adata.obs['Class'].isin(grouping[group])].copy()
    adata_sub.X = adata_sub.layers["UMIs"]
    ## Find HVGs
    sc.pp.highly_variable_genes(adata_sub, flavor="seurat_v3", n_top_genes=2000)
    ## Setup and train scVI model
    scvi.model.SCVI.setup_anndata(
        adata_sub,
        layer="UMIs",
        batch_key="donor_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
        categorical_covariate_keys=["donor_id"],
    )
    model = scvi.model.SCVI(adata_sub, 
                            dispersion="gene-batch", 
                            n_hidden=128, 
                            n_latent=64, 
                            n_layers=3)
    model.train(max_epochs=200)
    ## Get latent representation
    adata_sub.obsm["X_scVI"] = model.get_latent_representation()
    ## Compute UMAP
    sc.pp.neighbors(adata_sub, use_rep="X_scVI")
    sc.tl.umap(adata_sub)
    ## Plot UMAP
    colors = adata_sub.obs["Group"].map(color_map).values.astype(str)
    colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values
    ## Plot UMAP colored by cell type
    plt.figure(figsize=(12,12))
    plt.scatter(adata_sub.obsm["X_umap"][:,0], 
                adata_sub.obsm["X_umap"][:,1], 
                c=colors,
                s=0.1, 
                alpha=0.8)
    ## Add color legend
    plt.title("UMAP of HMBA BG Neurons - Consensus Annotation", fontsize=16)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "class_umaps", f"HMBA_consensus_{group}_UMAP.png"), dpi=900)
    ## Save the subsetted AnnData with UMAP
    outpath = os.path.join(work_dir, "subsets", f"HMBA_basalganglia_{group}_xspecies_UMAP.h5ad")
    adata_sub.write_h5ad(outpath)
    print(f"Saved UMAP for {group} to {outpath}")