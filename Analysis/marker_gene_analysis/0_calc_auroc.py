import os
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata


def compute_aurocs_markers(
    adata, 
    label_key="cell_type", 
    min_cells=10, 
    min_frac=0.2, 
    expr_thresh=1.0,
    auroc_thresh=0.9
):
    """
    Compute gene AUROCs with robust marker assignment.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object (cells × genes)
    label_key : str
        Column in adata.obs with cell type labels
    min_cells : int
        Minimum number of cells per cell type to include
    min_frac : float
        Minimum fraction of expressing cells to consider a gene as marker
    auroc_thresh : float
        AUROC threshold to call a gene a marker
    
    Returns
    -------
    result_df : pd.DataFrame
        AUROC DataFrame with added columns:
        - top_celltype
        - is_marker
    """
    X = adata.X
    labels = np.asarray(adata.obs[label_key])
    n_cells, n_genes = X.shape
    # Filter cell types with enough cells
    cell_counts = pd.Series(labels).value_counts()
    valid_celltypes = cell_counts[cell_counts >= min_cells].index
    valid_mask = np.isin(labels, valid_celltypes)
    X = X[valid_mask, :]
    labels = labels[valid_mask]
    # --- Rank transform
    if sparse.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X
    ranks = np.apply_along_axis(rankdata, 0, X_dense, method="average")
    # --- Indicator matrix for one-vs-all
    celltypes, label_codes = np.unique(labels, return_inverse=True)
    Y = sparse.csr_matrix((np.ones(len(labels)), (np.arange(len(labels)), label_codes)),
                          shape=(len(labels), len(celltypes)))
    # --- Sum of positive ranks
    sum_pos_ranks = Y.T @ ranks  # types × genes
    n_pos = np.asarray(Y.sum(axis=0)).ravel()
    n_neg = len(labels) - n_pos
    numer = sum_pos_ranks / n_pos[:, None] - (n_pos[:, None] + 1)/2
    aurocs = (numer / n_neg[:, None]).T  # genes × types
    aurocs_df = pd.DataFrame(aurocs, index=adata.var_names, columns=celltypes)
    # --- Compute fraction of expressing cells per gene × celltype
    expr_mask = X_dense >= expr_thresh
    expr_frac = pd.DataFrame(
        (expr_mask.T @ Y.toarray()) / n_pos[None, :],
        index=adata.var_names,
        columns=celltypes
    )
    ## --- Assign is_marker if AUROC and expression fraction pass thresholds
    mask = (aurocs_df[celltypes] >= auroc_thresh) & (expr_frac >= min_frac)
    aurocs_df["is_marker"] = mask.any(axis=1)
    ## --- Assign top celltype for marker genes only
    aurocs_df["top_celltype"] = aurocs_df[celltypes].idxmax(axis=1)
    aurocs_df.loc[aurocs_df["is_marker"] == False, "top_celltype"] = None
    return aurocs_df

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/"

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir, exist_ok=True)

## -------------------------------------------
## Identify species files, could also be a list from glob
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "Human_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "Marmoset_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
}

## ----------------------------------------
## Load cross species data
cross_species = {}
for k in species_h5ads.keys():
    print(f"Loading {k} data...")
    cross_species[k] = ad.read_h5ad(species_h5ads[k])

##
for anno_term in ["Class", "Subclass", "Group"]:
    ## Skip if already calculated
    # if os.path.exists(os.path.join(analysis_dir, "tau_scores", f"{anno_term.lower()}_tau_scores.csv")):
    #     continue
    print(f"Calculating AUROC for {anno_term}")
    ## Store AUROCs per species
    AUROCs = {}
    for species in species_h5ads:
        ## Load species data
        print(f"Processing {species}")
        ##
        AUROCs[species] = compute_aurocs_markers(cross_species[species], label_key=anno_term)
    ## Merge results using gene as the key
    aurocResults = pd.concat(AUROCs, axis=1)
    ## Fill NaNs with 0
    aurocResults = aurocResults.fillna(0)
    ## Add stats
    aurocResults['xspecies_mean'] = aurocResults.iloc[:,0:3].mean(axis=1)
    aurocResults['xspecies_median'] = aurocResults.iloc[:,0:3].median(axis=1)
    aurocResults['xspecies_min'] = aurocResults.iloc[:,0:3].min(axis=1)
    aurocResults['xspecies_max'] = aurocResults.iloc[:,0:3].max(axis=1)
    ## Arrange by decreasing xspecies_median
    aurocResults = aurocResults.sort_values(by='xspecies_min', ascending=False)
    # Step 5: Save to CSV
    aurocResults.to_csv(os.path.join(analysis_dir, f"{anno_term.lower()}_auroc_scores.csv"), index=True)






















# # --- logFC per gene × celltype (one-vs-all)
# logfc = pd.DataFrame(index=adata.var_names, columns=celltypes, dtype=float)
# for i, ct in enumerate(celltypes):
#     cells_in = (labels == ct)
#     cells_out = ~cells_in
#     mean_in = X_dense[cells_in, :].mean(axis=0)
#     mean_out = X_dense[cells_out, :].mean(axis=0)
#     logfc[ct] = np.log2(mean_in + 1) - np.log2(mean_out + 1)