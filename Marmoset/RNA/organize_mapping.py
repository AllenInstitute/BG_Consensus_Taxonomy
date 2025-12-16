import anndata as ad
import pandas as pd
import sciduck as sd

##
mapping_results = ad.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Marmoset/BasalGanglia/Marmoset_basalganglia_anno_latest_MapMyCells.h5ad", backed="r")

## Lets gather columns to put into annotation table
anno = mapping_results.obs.columns[mapping_results.obs.columns.str.contains("|".join(["AIT117_.*_label", "AIT193_.*_label", "Siletti_.*_label", "ABCmouse_.*_label", "HMBA_WB_Human_.*_label"]), regex=True)].tolist()
numeric_anno = mapping_results.obs.columns[mapping_results.obs.columns.str.contains("|".join([".*_avg_correlation", ".*_bootstrapping_probability"]), regex=True)].tolist()

##
anno_table = sd.anno.build_annotation_table(mapping_results, 
                                    group_by="cluster_id", 
                                    categorical_annotations=anno + ["donor_name"],
                                    numeric_annotations= numeric_anno, 
                                    min_percent=0.2, 
                                    annotation_alerts={"donor_name": 0.90})
anno_table.to_csv(os.path.join("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Marmoset/BasalGanglia/", f"Marmoset_consensus_anno_table.csv"))