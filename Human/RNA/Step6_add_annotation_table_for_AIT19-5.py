import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import sys
import sciduck as sd

##
species = "Human"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"
os.chdir(work_dir)


adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-5.h5ad")) 

## Lets gather columns to put into annotation table
anno = adata.obs.columns[adata.obs.columns.str.contains("|".join(["AIT117_.*_label", "AIT193_.*_label", "Silettisub_.*_label", "ABCmouse_.*_label", "HMBA_WB_Human_.*_label"]), regex=True)].tolist()
numeric_anno = adata.obs.columns[adata.obs.columns.str.contains("|".join([".*_avg_correlation", ".*_bootstrapping_probability"]), regex=True)].tolist()

## Previous work
anno_prev = ['Neighborhood_AIT19.3', 'Class_AIT19.3', 'Subclass_AIT19.3', 'Group_AIT19.3', 'cluster_AIT19.3']

## Additional metadata
meta = ["donor_name", "barcoded_cell_sample_label", "sex", "age", "anatomical_region", "ROI", "merged_ROIs"]

##
anno_table = sd.anno.build_annotation_table(adata, 
                                    group_by="cluster", 
                                    categorical_annotations=anno_prev + anno + meta,
                                    numeric_annotations=["doublet_score", "n_genes_by_counts", "umi.counts", "percent_keeper", "tso_percent"] + numeric_anno, 
                                    min_percent=0.05, 
                                    annotation_alerts={
                                        "donor_name": 0.90,
                                        "barcoded_cell_sample_label": 0.90,
                                        "sex": 0.85,
                                        "merged_ROIs": 0.9,
                                        "anatomical_region": 0.8    
                                    })
anno_table.to_csv(os.path.join(work_dir, f"{species}_{pipeline}_AIT19-5_consensus_anno_table.csv"))
