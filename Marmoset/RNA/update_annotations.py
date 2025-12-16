import anndata as ad
import pandas as pd
import os

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import cell_type_mapper
from cell_type_mapper.utils.anndata_utils import (read_df_from_h5ad)
from cell_type_mapper.taxonomy.utils import (get_taxonomy_tree, validate_taxonomy_tree)

##
species = "Marmoset"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables/"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies"

## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_anno_latest.h5ad"), backed="r")

## Load annotation sheet
species_anno = pd.read_excel(os.path.join(anno_dir, f"{species}_BG_annot_table.xlsx"))

##
adata.obs["Group_old"] = adata.obs["Group"].copy()

## Update annotations
for anno in ["Neighborhood", "Class", "Subclass", "Group"]:
    if anno not in adata.obs.columns:
        adata.obs[anno] = "Unknown"
    adata.obs[anno] = adata.obs[anno].astype(str)
    adata.obs[anno] = adata.obs.loc[:,"cluster_id"].replace(dict(zip(species_anno["cluster_id"], species_anno[anno].astype(str))))
    adata.obs[anno] = adata.obs[anno].astype("category")

##
taxonomy_tree = get_taxonomy_tree(
    obs_records=adata.obs.to_dict(orient='records'),
    column_hierarchy=["Neighborhood", "Class", "Subclass", "Group", "cluster_id"])

##
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_anno_041425.h5ad"))
adata.write_h5ad(os.path.join(cxgdir, f"{species}_basalganglia_anno_latest.h5ad"))

