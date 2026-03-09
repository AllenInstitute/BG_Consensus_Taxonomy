import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import os, sys, glob, re

## KNN Classifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript"
analysis_dir = os.path.join(base_dir, "Analysis", "species_specific_markers")
marker_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Accuracy_analysis"


## Ensure analysis dir exists
if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)

## Load gene sets
tf_database = pd.read_csv("/home/nelson.johansen/scripts/EvoGen/Projects/M1Evo/data/metadata/TFs_lambert_pmid29425488_1.01.txt", sep="\t")
tf_database = tf_database["hgnc_symbol"]

## Tau results 
marker_files = glob.glob(
    os.path.join(base_dir, "Analysis", "Expressologs", "tau_scores", "group_tau_scores_mean_donor_meta_analysis_all.csv")
)

markers = pd.read_csv(marker_files[0], index_col=0)

## Species column names
species_cols = ["human", "macaque", "mouse"]

## Boolean masks
high_mask = markers[species_cols] >= 0.95
low_mask  = markers[species_cols] < 0.05

## Find genes that satisfy: exactly one species high, all others low
exclusive_mask = (high_mask.sum(axis=1) == 1) & (low_mask.sum(axis=1) == (len(species_cols)-1))
exclusive_markers = markers[exclusive_mask]

## Add which species a gene is exclusive to
exclusive_markers = exclusive_markers.assign(
    marker_species=high_mask.idxmax(axis=1)  ## species with the high tau
)

exclusive_markers.marker_species.value_counts()

for species in species_cols:
    species_markers = exclusive_markers.loc[exclusive_markers.marker_species == species]
    ##
    with open(os.path.join(marker_dir, f"{species}_markers_092325.tsv"), "w") as f:
        for item in species_markers.index.tolist():
            f.write(str(item) + "\n")

## --------------------------------------------
## Plot the markers per species
## --------------------------------------------

##### -- Load Data -- ######
anno_term = "Group"
for species in species_cols:
    mean_expr = ad.read_h5ad(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Expressologs/{anno_term.lower()}/{species.lower()}_{anno_term.lower()}_mean_expr_all.h5ad")
    ## Subset to species specific markers
    species_markers = exclusive_markers.loc[exclusive_markers.marker_species == species]
    species_marker_set = [str(g) for g in species_markers.index.tolist() if g in mean_expr.var_names]
    ##
    adata_consensus_markers = mean_expr[:,species_marker_set]
    ## Merge species_markers info into adata.var
    adata_consensus_markers.var = adata_consensus_markers.var.merge(
        species_markers,
        left_index=True,
        right_index=True,
        how="left"
    )
    ## Save for plotting
    adata_consensus_markers.write_h5ad(os.path.join(analysis_dir, f"{species.lower()}_exclusive_tau_markers.h5ad"))

## --------------------------------------------
## Quick check on cell type accuracy with species specific markers
## --------------------------------------------

# ## Check human Group mapping accuracy
# adata = ad.read_h5ad(os.path.join(adata_dir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print.h5ad"))

# ## Unique markers to species
# marker_set = exclusive_markers.loc[exclusive_markers["marker_species"] == "Macaque"].index.tolist()
# marker_set = [g for g in marker_set if g in adata.var_names]

# ## Train KNN Classifier
# X = adata[:,marker_set].X 
# y = adata.obs["Group"]

# # KNN classifier
# knn = KNeighborsClassifier(n_neighbors=10, metric="euclidean")
# knn.fit(X, y)

# # Predictions
# y_pred = knn.predict(X)

# report_dict = classification_report(y, y_pred, output_dict=True)
# report = pd.DataFrame(report_dict).transpose()
