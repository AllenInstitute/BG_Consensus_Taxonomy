import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import pickle
import os

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
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9.h5ad"))

## Load in annotation transfer from MapMyCells
annotations = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_AIBS_BICAN_MapMyCells.h5ad"))
anno_columns = [anno for anno in annotations.obs.columns if "MapMyCells" in anno]
anno_columns = list(set(anno_columns) - set(adata.obs.columns))

## Transfer annotations
adata.obs = adata.obs.merge(annotations.obs[anno_columns], left_index=True, right_index=True)

## Load in annotation transfer from MapMyCells
annotations = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_AIBS_BICAN_MapMyCells_Flat.h5ad"))
anno_columns = [anno for anno in annotations.obs.columns if "MapMyCells" in anno]
anno_columns = list(set(anno_columns) - set(adata.obs.columns))

## Transfer annotations
adata.obs = adata.obs.merge(annotations.obs[anno_columns], left_index=True, right_index=True)

##
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno.h5ad")) 
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno.h5ad"))

