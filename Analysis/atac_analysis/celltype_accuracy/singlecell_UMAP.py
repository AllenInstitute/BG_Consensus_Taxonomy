import os
import sys
import glob
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import sparse
import scanpy as sc

## -------------------------  
##
## -------------------------
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
## Load data
## -------------------------
species = "human"
umap_df = pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/human/ATAC/h5ad/concatenated_peakmat_UMAP_harmony.csv', index_col=0)
umap_df["Group"] = umap_df["Group"].replace(" ", "_", regex=True)

## -------------------------
colors = umap_df["Group"].map(color_map).values.astype(str)
colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values
colors[colors == 'nan'] = "lightgrey"  # Handle NaN values

## Plot UMAP colored by cell type
plt.figure(figsize=(12,12))
plt.scatter(umap_df["UMAP1"], 
            umap_df["UMAP2"], 
            c=colors.tolist(),
            s=0.005, 
            alpha=0.8)
## Add color legend
plt.title("UMAP of variable Human peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "crested", "figures", f"{species}_snapatac2_umap_250000_features_harmony.png"), dpi=900)
