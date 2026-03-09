from __future__ import annotations

import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
import h5py
import os
import csv
import re
from scipy.sparse import csr_matrix

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
out_dir = os.path.join(analysis_dir, "atac", "specificity")

## -------------------------------------------
## Load Macaque chromosome alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## From CRESted
def calc_gini(targets: np.ndarray) -> np.ndarray:
    """
    Return Gini scores for the given targets.

    This function calculates the Gini coefficient for each row in the targets array
    and assigns the score to the maximum value's position in the corresponding row
    of the Gini scores array.

    Parameters
    ----------
    targets
        A 2D numpy array where each row represents a set of target values.

    Returns
    -------
    gini scores
        A 2D numpy array with the same shape as `targets` containing Gini scores,
        where each score is assigned to the position of the maximum value in each row.
    """
    def _gini(array: np.ndarray) -> float:
        """Calculate the Gini coefficient of a numpy array."""
        array = (
            array.flatten().clip(0, None) + 0.0000001
        )  # Ensure non-negative values and avoid zero
        array = np.sort(array)
        index = np.arange(1, array.size + 1)
        return (np.sum((2 * index - array.size - 1) * array)) / (
            array.size * np.sum(array)
        )
    ##
    gini_scores = np.zeros_like(targets)
    ##
    for region_idx in range(targets.shape[0]):
        region_scores = targets[region_idx]
        max_idx = np.argmax(region_scores)
        gini_scores[region_idx, max_idx] = _gini(region_scores)
    ##
    return gini_scores

species_h5ads = {
    "human" : os.path.join(base_dir, "human", "crested_adata", "human_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "crested_adata", "macaque_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "crested_adata", "marmoset_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
}

## -------------------------------------------
gini_scores_df = {}
for species, h5ad_path in species_h5ads.items():
    print(f"Processing {species} data...")
    adata = ad.read_h5ad(h5ad_path)
    ##
    if species == "human":
        adata = adata[~adata.obs_names.isin(["Alpha-MN", "Gamma-MN", "CSF-cNs"])].copy()  # Remove low number cell types
    ##
    if isinstance(adata.X, csr_matrix):
        target_matrix = (
            adata.X.toarray().T
        )  # Convert to dense and transpose to (regions, cell types)
    else:
        target_matrix = adata.X.T
    ##
    target_matrix = np.nan_to_num(target_matrix, nan=0.0)
    ##
    gini_scores = calc_gini(target_matrix)
    ##
    mean = np.mean(np.max(gini_scores, axis=1))
    std_dev = np.std(np.max(gini_scores, axis=1))
    gini_threshold = mean * std_dev
    ##
    adata.uns["gini_threshold"] = gini_threshold
    adata.layers["gini_scores"] = gini_scores.T
    adata.var["gini_scores"] = np.max(gini_scores, axis=1)
    adata.var["gini_scores"].to_csv(os.path.join(out_dir, species + "_gini_scores.csv"), index=True, header=True)
    ##
    gini_scores_df[species] = pd.DataFrame(adata.var["gini_scores"])
    gini_scores_df[species]["species"] = species
    gini_scores_df[species]["gini_threshold"] = gini_threshold
    ##
    adata.write_h5ad(os.path.join(out_dir, species + "_with_gini.h5ad"))

##
concatenated_df = pd.concat(gini_scores_df, axis=0)
concatenated_df.sort_values(by='gini_scores', ascending=False, inplace=True)

## Remove level 0 index
concatenated_df.reset_index(level=0, drop=True, inplace=True)
concatenated_df.to_csv(os.path.join(out_dir, "gini_scores_combined.csv"), index=True, header=True)