import numpy as np 
import keras
import crested
import anndata as ad
import os
import tensorflow as tf
import pysam
import pickle
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  ## non-GUI backend that supports tostring_rgb()
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

## -------------------------
## Helpful paths
## -------------------------
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
## Model related
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
rnadata_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/s3/AIT/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
## Data
## -------------------------
species = "human"
adata = ad.read_h5ad(os.path.join(data_dir, f"{species}/crested_adata/{species}_basalganglia_hmba_pre-print_finetune_crested.h5ad"))
adata_rna = os.path.join(rnadata_dir, f"{species.capitalize()}_HMBA_basalganglia_AIT_pre-print.h5ad")

## -------------------------
## Plot pattern matrix
## -------------------------

# crested.tl.contribution_scores_specific(
#     input=adata_filtered,
#     target_idx=None,  # We calculate for all classes
#     model=model,
#     method='integrated_grad',
#     output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/patterns/modisco_results_ft_500/",
#     verbose=True
# )

# import crested
# meme_db, motif_to_tf_file = crested.get_motif_db()

# import os
# os.environ["PATH"] = "/home/niklas.kempynck/meme/bin:" + os.environ["PATH"]

# # run tfmodisco on the contribution scores
# crested.tl.modisco.tfmodisco(
#     window=1000,
#     output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/patterns/modisco_results_ft_500/",
#     contrib_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/patterns/modisco_results_ft_500/",
#     #report=True,  # Optional, will match patterns to motif MEME database
#     #meme_db=meme_db,  # File to MEME database
#     max_seqlets=20000,
# )

## First we obtain the resulting modisco files per class
matched_files = crested.tl.modisco.match_h5_files_to_classes(
    contribution_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/patterns/modisco_results_ft_500/", classes=list(adata.obs_names)
)

sim_matrix, pattern_ids, pattern_dict = crested.tl.modisco.calculate_tomtom_similarity_per_pattern(
    matched_files = matched_files,
    trim_ic_threshold=0.025,
    verbose=True
)

groups = []
for id in pattern_ids:
    ct = '_'.join(id.split('_')[:-3])
    groups.append(ct)

unique_cats = pd.unique(groups)
group_colors = {cat: to_hex(get_cmap("tab20", len(unique_cats))(i)) for i, cat in enumerate(unique_cats)}

##
with plt.rc_context({'figure.figsize': (24,24), 'figure.dpi': 900}):
    crested.pl.patterns.clustermap_tomtom_similarities(
        sim_matrix=sim_matrix,
        ids=pattern_ids,
        pattern_dict=pattern_dict,
        group_info = [(groups, group_colors)], # Grouping labels
        min_seqlets=200, # Add a minimum amount of seqlets to take the most relevant patterns
        figsize=(24,24),
        save_path=os.path.join(analysis_dir, f"crested/patterns/figures/{species}_finetune_tomtom_pattern_similarity_clustermap.pdf")
    )

## -------------------------
## Create cluster map with pwm logos
## -------------------------

## Load in patterns
with open("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/patterns/patterns_updated.pkl", "rb") as f:   
    all_patterns = pickle.load(f)

## Create pattern matrix
pattern_matrix = crested.tl.modisco.create_pattern_matrix(
    classes=list(adata.obs_names),
    all_patterns=all_patterns,
    normalize=False,
    pattern_parameter="seqlet_count_log",
)
pattern_matrix.shape

##
with plt.rc_context({'figure.figsize': (30,10), 'figure.dpi': 600}):
    crested.pl.patterns.clustermap_with_pwm_logos(
        pattern_matrix,
        list(adata.obs_names),
        pattern_dict=all_patterns,
        figsize=(30, 10),
        grid=True,
        dendrogram_ratio=(0.03, 0.08),
        importance_threshold=4.5,
        logo_height_fraction=0.35,
        logo_y_padding=0.1,
        save_path=os.path.join(analysis_dir, f"crested/patterns/figures/{species}_finetune_model_clustermap_with_pwm_logos.pdf")
    )

## -------------------------
## Create cluster map with pwm logos
## -------------------------

# mean_expression_df = crested.tl.modisco.calculate_mean_expression_per_cell_type(
#     adata_rna, "Group", cpm_normalize=False
# )
# mean_expression_df.to_csv(os.path.join(analysis_dir, f"crested/patterns/data/{species}_mean_expression_per_cell_type.csv"))
mean_expression_df = pd.read_csv(os.path.join(analysis_dir, f"crested/patterns/data/{species}_mean_expression_per_cell_type.csv"), index_col=0)

##
classes = cluster_meta["Group"][cluster_meta["Group"].isin(adata.obs_names)].values

##
mean_expression_df.index = [s.replace(" ", "_") for s in mean_expression_df.index.values]
mean_expression_df = mean_expression_df.loc[classes]

##
contribution_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/{species}/patterns/modisco_results_ft_500/"
html_paths = crested.tl.modisco.generate_html_paths(
    all_patterns, classes, contribution_dir
)

## q_val threshold to only select significant matches
pattern_match_dict = crested.tl.modisco.find_pattern_matches(
    all_patterns, html_paths, q_val_thr=0.05
)

## Get motif to tf mapping
meme_db, motif_to_tf_file = crested.get_motif_db()
motif_to_tf_df = crested.tl.modisco.read_motif_to_tf_file(
   motif_to_tf_file
)

cols = [
    "Mouse_Direct_annot",
    "Mouse_Orthology_annot",
    "Cluster_Mouse_Direct_annot",
    "Cluster_Mouse_Orthology_annot",
]
cols = [
    "Human_Direct_annot",
    "Human_Orthology_annot",
    "Cluster_Human_Direct_annot",
    "Cluster_Human_Orthology_annot",
]

##
pattern_tf_dict, all_tfs = crested.tl.modisco.create_pattern_tf_dict(
    pattern_match_dict, motif_to_tf_df, all_patterns, cols
)
tf_ct_matrix, tf_pattern_annots = crested.tl.modisco.create_tf_ct_matrix(
    pattern_tf_dict,
    all_patterns,
    mean_expression_df,
    classes,
    log_transform=False,
    normalize_pattern_importances=False,
    normalize_gex=True,
    min_tf_gex=0.95,
    importance_threshold=4,
    pattern_parameter="seqlet_count_log",
    filter_correlation=True,
    verbose=True,
    zscore_threshold=1.5,
    correlation_threshold=0.5,
)

##
with plt.rc_context({'figure.figsize': (25,10), 'figure.dpi': 900}):
    crested.pl.patterns.clustermap_tf_motif(
        tf_ct_matrix,
        heatmap_dim="contrib",
        dot_dim="gex",
        class_labels=classes,
        pattern_labels=tf_pattern_annots,
        fig_size=(25, 10),
        cluster_rows=False,
        cluster_columns=True,
        save_path=os.path.join(analysis_dir, f"crested/patterns/figures/{species}_finetune_tf_ct_clustermap.pdf")
    )