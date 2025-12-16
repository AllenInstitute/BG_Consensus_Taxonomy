import anndata as ad
import pandas as pd
import os

import cell_type_mapper
from cell_type_mapper.utils.anndata_utils import (read_df_from_h5ad)
from cell_type_mapper.taxonomy.utils import (get_taxonomy_tree, validate_taxonomy_tree)

##
species = "Macaque"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies"

## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno_latest.h5ad"))

## Load annotation sheet
macaque_anno = pd.read_excel(os.path.join(anno_dir, f"{species}_AIBS_AIT11-9_consensus_anno_table.xlsx"))
macaque_anno.cluster_id = "Macaque-" + macaque_anno.cluster_id.astype(str)

## Update Macaque
for anno in ["Neighborhood", "Class", "Subclass", "Group"]:
    if anno not in adata.obs.columns:
        adata.obs[anno] = "Unknown"
    adata.obs[anno] = adata.obs[anno].astype(str)
    adata.obs[anno] = adata.obs.loc[:,"cluster_id"].replace(dict(zip(macaque_anno["cluster_id"], macaque_anno[anno].astype(str))))
    adata.obs[anno] = adata.obs[anno].astype("category")

##
taxonomy_tree = get_taxonomy_tree(
    obs_records=adata.obs.to_dict(orient='records'),
    column_hierarchy=["Neighborhood", "Class", "Subclass", "Group", "cluster_id"])

##
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno_042125.h5ad"))
adata.write_h5ad(os.path.join(cxgdir, f"{species}_basalganglia_HMBA_AIT11-9_anno_latest.h5ad"))

##
adata = adata[:,adata.var_names[adata.var.highly_variable]]
adata.layers["UMIs"] = adata.layers["UMIs"].astype("float32") ## For R..., dgR matrix must be double precision
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno_042125_AIT.h5ad"))


