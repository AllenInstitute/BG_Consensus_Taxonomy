import pandas as pd
import numpy as np 
from scipy.sparse import issparse
import scanpy as sc

work_dir = "./" # please change it to "/your/desired/output/path"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/BICAN-releases/final"
species = "Human"

adata = sc.read_h5ad(f"{data_dir}/{species}_HMBA_basalganglia_AIT_pre-print.h5ad")


# Set expression threshold (log2(CPM) > 3 or log1p(CPM) > 2.197 depending on scale)
EXPRESSION_THRESHOLD = 1  # 2.197 seems too strict in human BG; 1 still ?!
CLUSTER_KEY = 'Cluster'  # Key in .obs that defines cell clusters


# Define neurotransmitter marker logic (each value is a list of OR-groups; multiple groups imply AND logic)
NEUROTRANSMITTER_MARKERS = {
    "Glut": [["SLC17A6", "SLC17A7", "SLC17A8"]],
    "GABA": [["SLC32A1", "SLC18A2"], ["GAD1", "GAD2", "ALDH1A1"]],
    "Glyc": [["SLC6A5"]],
    "Chol": [["SLC18A3"], ["CHAT"]],
    "Dopa": [["SLC6A3", "SLC18A2"], ["TH", "DDC"]],
    "Sero": [["SLC6A4", "SLC18A2"], ["TPH2"], ["DDC"]],
    "Nora": [["SLC6A2", "SLC18A2"], ["DBH"]],
    "Hist": [["SLC18A2"], ["HDC"]],
}


# Flatten all genes used in the marker logic
all_marker_genes = sorted({gene for groups in NEUROTRANSMITTER_MARKERS.values() for group in groups for gene in group})
genes_available = [g for g in all_marker_genes if g in adata.var_names]


# Extract expression matrix for marker genes only
gene_indices = [adata.var_names.get_loc(g) for g in genes_available]
expr_df = pd.DataFrame(
    adata.X[:, gene_indices].toarray(),
    columns=genes_available,
    index=adata.obs_names
)
expr_df[CLUSTER_KEY] = adata.obs[CLUSTER_KEY].values


# Compute median expression per cluster
cluster_medians = expr_df.groupby(CLUSTER_KEY).median()


# Define logic function & apply logic to each cluster, return list of matching neurotransmitters
def check_nt_identity(row, marker_logic, threshold=3.0):
    results = []
    for nt, rules in marker_logic.items():
        satisfied = True
        for or_group in rules:
            if not any((gene in row.index and row[gene] > threshold) for gene in or_group):
                satisfied = False
                break
        if satisfied:
            results.append(nt)
    return results

nt_identity = cluster_medians.apply(
    lambda row: check_nt_identity(row, NEUROTRANSMITTER_MARKERS, threshold=EXPRESSION_THRESHOLD),
    axis=1
)

# Convert list to boolean matrix
nt_identity_bool_df = pd.DataFrame(False, index=nt_identity.index, columns=NEUROTRANSMITTER_MARKERS.keys())

for cluster, matched_types in nt_identity.items():
    for nt in matched_types:
        nt_identity_bool_df.loc[cluster, nt] = True

nt_identity_bool_df.insert(0, 'cluster', nt_identity_bool_df.index)

# csv output 
nt_identity_bool_df.to_csv(
    f"{work_dir}/{species}_BG_neurotransmitter_identity_by_cluster_th{EXPRESSION_THRESHOLD}.csv",
    index=False
)

print(nt_identity_bool_df)
