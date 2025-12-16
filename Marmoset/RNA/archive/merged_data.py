import scanpy as sc
import pandas as pd
import glob
import os

##
os.chdir('/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Marmoset/raw_data/marm_hmba_h5ad')

## Gather all .h5ad file names
marmoset_files = glob.glob("*.h5ad")

## Read in all .h5ad files
marmoset_adata = []
for file in marmoset_files:
    adata = sc.read_h5ad(file)
    marmoset_adata.append(adata)

## Concatenate all .h5ad files
marmoset_adata = sc.AnnData.concatenate(*marmoset_adata)

## Metadata
meta = pd.read_csv("datalog_metadata.csv")

## Add metadata
merged_metadata = marmoset_adata.obs.merge(meta, left_on="sample_name", right_on="sample_name")
merged_metadata.index = marmoset_adata.obs.index
marmoset_adata.obs = merged_metadata

## Save
marmoset_adata.write("marmoset_HMBA_BG_merged.h5ad")


