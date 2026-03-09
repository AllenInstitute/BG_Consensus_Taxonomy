import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, glob, re

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"

## Load gene sets
tf_database = pd.read_csv("/home/nelson.johansen/scripts/EvoGen/Projects/M1Evo/data/metadata/TFs_lambert_pmid29425488_1.01.txt", sep="\t")
tf_database = tf_database["hgnc_symbol"]

## ----------------------------------------
## Tau results
## ----------------------------------------
marker_files = glob.glob(
    os.path.join(analysis_dir, "marker_gene_analysis", "tau_scores/*_tau_scores_mean_donor_meta_analysis_all.csv")
)

## Gather up most conserved markers across annotation levels
marker_collection = []
for file in marker_files:
    ## 
    markers = pd.read_csv(file, index_col=0)
    ## Apply threshold to group and then higher level annotations
    if re.search("group", file):
        filtered = markers[markers["xspecies_min"] > 0.70]
    else:
        filtered = markers[markers["xspecies_min"] > 0.9]
    ##
    marker_collection.append(set(filtered.index))

## Union of all sets
marker_set = set().union(*marker_collection)

## -------------------------------------------
## Expressolog / AUROC results
## -------------------------------------------

##
expressologs = pd.read_csv(os.path.join(analysis_dir, "expressologs", "group", "group_expressologs_with_metrics.csv"), index_col=0)

## Gather conserved markers
conserved_tf_markers = (
    expressologs
    .loc[expressologs.index.isin(marker_set)]
    .loc[lambda df: df.index.isin(tf_database)]
    .loc[lambda df: df["auroc"] > 0.70]
    .loc[lambda df: ~df.index.duplicated(keep='first')]
    .index
    .tolist()
)

conserved_markers = (
    expressologs
    .loc[expressologs.index.isin(marker_set)]
    .loc[lambda df: df["auroc"] > 0.70]
    .loc[lambda df: ~df.index.duplicated(keep='first')]
    .index
    .tolist()
)

with open(os.path.join(analysis_dir, "marker_gene_analysis", "conserved_markers_092325.tsv"), "w") as f:
    for item in conserved_markers:
        f.write(str(item) + "\n")

with open(os.path.join(analysis_dir, "marker_gene_analysis", "conserved_tf_markers_092325.tsv"), "w") as f:
    for item in conserved_tf_markers:
        f.write(str(item) + "\n")

## -------------------------------------------
## Marker data summary
## -------------------------------------------
adata_consensus = ad.read_h5ad(os.path.join(base_dir, "xspecies", "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"))
adata_consensus.X = adata_consensus.layers["UMIs"]

sc.pp.normalize_total(adata_consensus, target_sum=1e6)
sc.pp.log1p(adata_consensus)

## Merge into a summary by Group
adata_consensus_markers = adata_consensus[:,conserved_markers]
adata_consensus_markers = sc.get.aggregate(adata_consensus_markers, by="Group", func=["mean"])
adata_consensus_markers.write_h5ad(os.path.join(analysis_dir, "marker_gene_analysis", "AIBS_BG_consensus_taxonomy_harmonized_AIT-pre-print_conserved_markers_10925.h5ad"))

## Merge into a summary by Group
adata_consensus_tf = adata_consensus[:,adata_consensus.var_names.isin(conserved_tf_markers)]
adata_consensus_tf = sc.get.aggregate(adata_consensus_tf, by="Group", func=["mean"])
adata_consensus_tf.write_h5ad(os.path.join(analysis_dir, "marker_gene_analysis", "AIBS_BG_consensus_taxonomy_harmonized_AIT-pre-print_tf_markers_10925.h5ad"))
