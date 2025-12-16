import anndata as ad
import sciduck as sd
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

## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno.h5ad"), backed="r")

## Lets gather columns to put into annotation table
anno = adata.obs.columns[adata.obs.columns.str.contains("|".join(["AIT117_.*_label", "AIT193_.*_label", "Siletti_.*_label", "ABCmouse_.*_label", "HMBA_WB_Human_.*_label"]), regex=True)].tolist()
numeric_anno = adata.obs.columns[adata.obs.columns.str.contains("|".join([".*_avg_correlation", ".*_bootstrapping_probability"]), regex=True)].tolist()

## Previous work
anno_prev = ["AIT117_Neighborhood", "AIT117_Class", "AIT117_Subclass", "AIT117_Group", "AIT117_cluster"]

## Additional metadata
meta = ["donor_name", "anatomical_region", "sex", "age", "barcoded_cell_sample_label"]

##
anno_table = sd.anno.build_annotation_table(adata, 
                                    group_by="cluster", 
                                    categorical_annotations=anno_prev + anno + meta,
                                    numeric_annotations=["doublet_score", "total_genes", "total_counts", "percent_keeper"] + numeric_anno, 
                                    min_percent=0.05, 
                                    annotation_alerts={"donor_name": 0.90})
anno_table.to_csv(os.path.join(work_dir, f"{species}_{pipeline}_consensus_anno_table.csv"))