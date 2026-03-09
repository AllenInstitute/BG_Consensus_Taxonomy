import os
import anndata as ad
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
import pertpy as pt
import scanpy as sc

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
adata = ad.read_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_neuron_alignment_labeltransfer_allgenes.h5ad"))

## -------------------------
adata_mouse = adata[adata.obs["organism"] == "Mouse", :].copy()

anno_term = "Group_propagated"
species = "mouse"
## Compute mean expression per Subclass
try:
    mean_expr = sc.get.aggregate(adata_mouse, by=anno_term, func=["mean"])
    mean_expr.write_h5ad(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/mean_expr/group/{species.lower()}_{anno_term.lower()}_mean_expr_orthologs.h5ad")
except:
    print(f"Error in {species}")
