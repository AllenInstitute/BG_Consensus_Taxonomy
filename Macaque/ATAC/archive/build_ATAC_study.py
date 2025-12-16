import scanpy as sc
import anndata as ad
import os

##
os.chdir("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque")

##
adata = sc.read_h5ad("HMBA_Macaque_BG_082024.h5ad", backed="r")

## Pull meta.data relevant to ATAC fragment files
study_set = adata.obs.drop_duplicates(subset=["ar_id"]).loc[:, ["species", "ar_id", "ar_directory"]]
study_set.to_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC/ATAC_HMBA_Macaque_BG_study.csv", index=False)