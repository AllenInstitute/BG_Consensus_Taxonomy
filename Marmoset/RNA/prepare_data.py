import anndata as ad
import pandas as pd
import os

##
species = "Marmoset"

## Helpful locations which are assumed to already exist
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies"

## Read in marmoset data
adata = ad.read_h5ad(os.path.join(work_dir, "rna_bg_marm.h5ad"))

##
adata.X = adata.layers["UMIs"].copy()

## Gather consistent metadata
adata.obs["cell_barcode"] = [cell_id.split('_')[-1].split('-')[0] for cell_id in adata.obs_names.tolist()]  # Extract cell barcode from index

##
adata.obs.rename(columns={"Cluster": "cluster_id", 
                            "Marm_Group_v3": "Group",
                            "barcoded_cell_sample_name": "barcoded_cell_sample_label",}, 
                            inplace=True)


adata.obs["cell_label"] = adata.obs["cell_barcode"].astype(str) + "-" + adata.obs["barcoded_cell_sample_label"].astype(str)
adata.obs["alignment_job_database"] = "Princeton"

## donor info
def extract_terms(elements):
    return [e.split("_")[2] if len(e.split("_")) > 2 else None for e in elements]

adata.obs["donor_id"] = "unknown"
adata.obs["donor_id"] = extract_terms(adata.obs["orig.ident"].tolist())
adata.obs["donor_id"] = adata.obs["donor_id"].astype('category')

## 
adata.obs.index = adata.obs["cell_label"].copy()

##
adata.raw = adata

## Write file
adata.write_h5ad(os.path.join(work_dir, "Marmoset_basalganglia_anno_latest.h5ad"))