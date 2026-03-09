import os
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp

def calc_tau(mat: pd.DataFrame):
    """
    Calculate Tau scores per gene from a gene x cell_type DataFrame.

    Parameters:
    mat (pd.DataFrame): DataFrame with genes as rows and cell types as columns.

    Returns:
    pd.Series: Tau score for each gene (index = gene names).
    """
    ## Get max expression per gene (row-wise)
    row_max = mat.max(axis=1)
    ## Avoid division by zero
    mat_norm = mat.div(row_max, axis=0)
    ## Compute Tau: sum(1 - norm) / (n_celltypes - 1)
    tau = (1 - mat_norm).sum(axis=1) / (mat.shape[1] - 1)
    ## Replace any NaNs with 0
    tau = tau.fillna(0)
    return tau

def get_cl_prop(adata, cluster_labels, min_cells = 10, expr_threshold=1, prop_threshold=0.5):
    """
    Compute cluster-level expression proportions for each gene.

    Parameters:
    adata (anndata.AnnData): Cell x gene matrix (assumed sparse).
    cluster_labels (array-like): Cluster label for each cell (same length as n_obs in adata).
    expr_threshold (float): Minimum expression value to consider a gene "expressed".
    prop_threshold (float): Proportion threshold to set gene expression to zero.

    Returns:
    pd.DataFrame: Gene x cluster matrix with expression proportions.
    """
    X = adata.X  # cell x gene matrix
    if not sp.issparse(X):
        raise ValueError("Expected a sparse matrix in adata.X")
    ##
    cluster_labels = np.asarray(cluster_labels)
    unique_clusters = np.unique(cluster_labels)
    # Create binary matrix: 1 if gene is expressed above threshold, 0 otherwise
    binarized = (X > expr_threshold).astype(np.int32)  # still sparse
    # Create cluster indicator matrix: n_cells x n_clusters
    from scipy.sparse import csr_matrix
    cluster_map = {c: i for i, c in enumerate(unique_clusters)}
    cl_indices = np.vectorize(cluster_map.get)(cluster_labels)
    n_clusters = len(unique_clusters)
    n_cells = adata.n_obs
    cluster_matrix = csr_matrix(
        (np.ones(n_cells), (np.arange(n_cells), cl_indices)),
        shape=(n_cells, n_clusters)
    )
    ## Compute gene x cluster counts using tcrossprod (binarized.T @ cluster_matrix)
    gene_cluster_counts = binarized.T.dot(cluster_matrix)  # shape: n_genes x n_clusters
    ## Normalize by number of cells in each cluster
    cluster_sizes = np.array(cluster_matrix.sum(axis=0)).ravel()  # shape: (n_clusters,)
    gene_cluster_prop = gene_cluster_counts.toarray() / cluster_sizes
    ## Set proportion to zero if gene is expressed in less than X% of cells for that cluster
    gene_cluster_prop[gene_cluster_prop < prop_threshold] = 0.0
    # max_props = gene_cluster_prop.max(axis=1)  # shape: n_genes,
    # mask = max_props < prop_threshold
    # gene_cluster_prop[mask, :] = 0.0
    ## Set clusters with min cells to zero
    gene_cluster_prop[:, cluster_sizes < min_cells] = 0.0
    ## Return as DataFrame with genes and cluster names
    return pd.DataFrame(
        gene_cluster_prop,
        index=adata.var_names,
        columns=unique_clusters
    )

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/"

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir, exist_ok=True)

if not os.path.exists(os.path.join(analysis_dir, "tau_scores")):
    os.makedirs(os.path.join(analysis_dir, "tau_scores"), exist_ok=True)

if not os.path.exists(os.path.join(analysis_dir, "tau_scores", "donor_level")):
    os.makedirs(os.path.join(analysis_dir, "tau_scores", "donor_level"), exist_ok=True)

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## 
celltype_set = cluster_meta.loc[cluster_meta.Class == "CN LGE GABA", "Group"].unique().tolist()

## Remove GPe_MEIS2-SOX6_GABA
celltype_set.remove("GPe_MEIS2-SOX6_GABA")

## -------------------------------------------
## Identify species files, could also be a list from glob
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "Human_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "Marmoset_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
}

## 
tau_threshold = 0.9
num_donor_sampling = 10

## ----------------------------------------
## Load cross species data
species_adatas = {}
for species in ["human", "macaque", "marmoset"]: 
    print(f"Loading {species} data...")
    adata = ad.read_h5ad(species_h5ads[species])
    adata.obs["Group"] = adata.obs["Group"].astype(str)
    adata.obs["Group"] = adata.obs["Group"].replace(" ", "_", regex=True)
    adata = adata[adata.obs["Group"].isin(celltype_set), :].copy()
    print(adata.shape)
    species_adatas[species] = adata

## Remove raw for each adata as we don't need it
for species in species_adatas:
    if species_adatas[species].raw is not None:
        species_adatas[species].raw = None



## Calculate Tau scores
for anno_term in ["Subclass", "Group"]:
    ## Skip if already calculated
    # if os.path.exists(os.path.join(analysis_dir, "tau_scores", f"{anno_term.lower()}_tau_scores.csv")):
    #     continue
    print(f"Calculating Tau scores for {anno_term}")
    tauScores_summary = {}; tauRecurrence_summary = {}; tauScores = {}
    for species in species_adatas:
        ## Load species data
        print(f"Processing {species}")
        tauScores[species] = pd.DataFrame(index=species_adatas[species].var_names)
        ## Number of donors to sample
        num_donor = max(2, int(np.ceil(len(species_adatas[species].obs['donor_id'].unique())/3)))
        for itr in range(num_donor_sampling):
            print(f"Iteration {itr+1} of {num_donor_sampling}")
            ## Gather sampled donors
            sampled_donors = np.random.choice(species_adatas[species].obs['donor_id'].unique(), size=num_donor, replace=False)
            donor_adata = species_adatas[species][species_adatas[species].obs['donor_id'].isin(sampled_donors)]
            ## Compute cluster-level proportions
            # propExpr = get_cl_prop(donor_adata, donor_adata.obs[anno_term])
            ## Means
            propExpr = sc.get.aggregate(donor_adata, by=anno_term, func=["mean"])
            ## Filter to genes with mean expression > 1 in at least one cluster
            propExpr = propExpr[:, (propExpr.layers["mean"] > 1).any(axis=0)]
            ## Convert to DataFrame with genes as rows, cell types as columns
            propExpr = pd.DataFrame(propExpr.layers["mean"].T, index=propExpr.var_names, columns=propExpr.obs[anno_term])
            ## Get tau scores
            tauScore = calc_tau(propExpr)
            tauScore = tauScore.sort_values(ascending=False)
            ## Save tau scores
            tauScores[species][f"{itr}_tau"] = tauScore
        ## Now compute mean and recurrence (# donors with tau > 7)
        tauScores[species]['mean_tau'] = tauScores[species].mean(axis=1)
        tauScores[species]['recurrence'] = (tauScores[species].iloc[:,0:-1] > tau_threshold).sum(axis=1)
        tauScores[species].to_csv(os.path.join(analysis_dir, "tau_scores", "donor_level", f"{species.lower()}_{anno_term.lower()}_tau_scores.csv"), index=True)
        ##
        tauScores_summary[species] = tauScores[species]['mean_tau'].copy()
        tauRecurrence_summary[species] = tauScores[species]['recurrence'].copy()
    ## Merge results using gene as the key
    tauResults = pd.concat(tauScores_summary, axis=1)
    tauResults.columns = tauResults.columns.get_level_values(0) ## Keep only species name
    ## Fill NaNs with 0
    tauResults = tauResults.fillna(0)
    ## Add stats
    tauResults['xspecies_mean'] = tauResults.iloc[:,0:3].mean(axis=1)
    tauResults['xspecies_median'] = tauResults.iloc[:,0:3].median(axis=1)
    tauResults['xspecies_min'] = tauResults.iloc[:,0:3].min(axis=1)
    tauResults['xspecies_max'] = tauResults.iloc[:,0:3].max(axis=1)
    ## Add in recurrence
    tauRecurrence = pd.concat(tauRecurrence_summary, axis=1)
    tauRecurrence.columns = tauRecurrence.columns.get_level_values(0) ## Keep only species name
    tauRecurrence = tauRecurrence.fillna(0)
    tauResults['xspecies_recurrence'] = tauRecurrence.mean(axis=1)
    ## Arrange by decreasing xspecies_median
    tauResults = tauResults.sort_values(by='xspecies_min', ascending=False)
    # Step 5: Save to CSV
    tauResults.to_csv(os.path.join(analysis_dir, "tau_scores", "MSN_vignette", f"{anno_term.lower()}_tau_scores_mean_donor_meta_analysis_all.csv"), index=True)
