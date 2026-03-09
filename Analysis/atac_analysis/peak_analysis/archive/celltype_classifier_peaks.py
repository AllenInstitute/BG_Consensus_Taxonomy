import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import f1_score, roc_auc_score, classification_report
from sklearn.preprocessing import label_binarize
from tqdm import tqdm

from __future__ import annotations

import anndata as ad
import os

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

import snapatac2 as snap

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
out_dir = os.path.join(analysis_dir, "atac", "specificity")

## -------------------------------------------
## Load Macaque chromosome alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## AnnoTable
cluster_meta = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## From CRESted
species_h5ads = {
    "human" : os.path.join(out_dir, "human_with_gini.h5ad"),
    "macaque" : os.path.join(out_dir, "macaque_with_gini.h5ad"),
    "marmoset": os.path.join(out_dir, "marmoset_with_gini.h5ad"),
}

def classify_cell_types(adata, group_by, n_permutations=100, random_state=42):
    """
    Quantify how accurately each cell type can be classified from a set of peaks.
    
    Parameters
    ----------
    X : pd.DataFrame or np.ndarray
        Cells x informative peaks matrix
    y : pd.Series or np.ndarray
        Cell type labels
    n_permutations : int
        Number of permutations for significance testing
    random_state : int
        Random seed for reproducibility
    
    Returns
    -------
    results_df : pd.DataFrame
        Per-cell-type metrics including F1-score, AUROC, and p-value
    """
    ## Prepare data
    rng = np.random.default_rng(random_state)
    X = adata.X
    y = adata.obs[group_by].values
    ## Fit classifier with cross-validation
    clf = LogisticRegression(multi_class='multinomial', solver='lbfgs', max_iter=1000)
    y_pred = cross_val_predict(clf, X, y, cv=5)
    ## Compute per-cell-type F1-score
    report = classification_report(y, y_pred, output_dict=True)
    f1_scores = {ct: report[ct]['f1-score'] for ct in np.unique(y)}
    ## Compute AUROC per cell type (one-vs-rest)
    y_bin = label_binarize(y, classes=np.unique(y))
    y_scores = cross_val_predict(clf, X, y, cv=5, method='predict_proba')
    auc_scores = {}
    for i, ct in enumerate(np.unique(y)):
        try:
            auc = roc_auc_score(y_bin[:, i], y_scores[:, i])
        except ValueError:
            auc = np.nan  # handle rare cases
        auc_scores[ct] = auc
    ## Permutation test for significance
    perm_f1 = {ct: [] for ct in np.unique(y)}
    for _ in tqdm(range(n_permutations), desc="Permutation test"):
        y_shuffled = rng.permutation(y)
        y_pred_perm = cross_val_predict(clf, X, y_shuffled, cv=5)
        report_perm = classification_report(y_shuffled, y_pred_perm, output_dict=True)
        for ct in np.unique(y):
            perm_f1[ct].append(report_perm[ct]['f1-score'])
    ## Compute p-values
    p_values = {}
    for ct in np.unique(y):
        perm_scores = np.array(perm_f1[ct])
        observed = f1_scores[ct]
        # one-sided test: how often permuted >= observed
        p = (np.sum(perm_scores >= observed) + 1) / (n_permutations + 1)
        p_values[ct] = p
    ## Combine results
    results_df = pd.DataFrame({
        'F1_score': f1_scores,
        'AUROC': auc_scores,
        'p_value': p_values
    })
    return results_df

## --------------------------------------------------
## Gather featurs for classifier
## --------------------------------------------------
h5ad_path = species_h5ads["macaque"]
print(f"Processing {species}...")

species_adata = ad.read_h5ad(h5ad_path)
del species_adata.obs["file_path"]

## Identify most accessible cell type for each peak
max_idx = np.argmax(species_adata.X, axis=0)
max_cell_types = species_adata.obs_names[max_idx]
species_adata.var["max_cell_type"] = max_cell_types.tolist()

## Gather top 100 peaks per cell type based on specificity (Gini)
top_peaks = []
for cell_type in species_adata.obs_names.unique():
    cell_type_peaks = species_adata.var[species_adata.var["max_cell_type"] == cell_type]
    top_cell_type_peaks = cell_type_peaks.nlargest(100, "gini_scores")
    top_peaks.append(top_cell_type_peaks.index.tolist())

##
top_peaks = [peak for sublist in top_peaks for peak in sublist]  # Flatten list
top_peaks = list(set(top_peaks))  # Unique peaks

##
adata_feature = species_adata[:,top_peaks].copy()

