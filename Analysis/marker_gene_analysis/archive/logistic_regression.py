import os, glob, re
import pandas as pd
import numpy as np
import anndata as ad

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

def species_celltype_markers_filtered(adatas, target_species, target_cell,
                                      top_k=200, expr_thresh=0.5):
    """
    Find species-specific cell type marker genes using Random Forest with optional filtering.
    
    Parameters
    ----------
    adatas : dict
        Keys = species, values = AnnData objects (cells x genes)
    target_species : str
        Species of interest
    target_cell : str
        Cell type of interest
    top_k : int
        Number of top genes to consider from classifier
    expr_thresh : float
        Maximum mean expression allowed in other species for same cell type (filter)
    
    Returns
    -------
    species_specific_genes : list
        Genes that are predictive and species-specific
    importances : pd.Series
        Global feature importances from Random Forest
    """
    # Step 1: concatenate data
    X_list, y_list, var_names = [], [], None
    for sp, ad in adatas.items():
        X_list.append(ad.layers["mean"])
        y_list.extend([f"{sp}_{ct}" for ct in ad.obs['Group']])
        if var_names is None:
            var_names = ad.var_names
    X = np.vstack(X_list)
    y = np.array(y_list)
    # Encode labels
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    target_label = f"{target_species}_{target_cell}"
    target_idx = np.where(le.classes_ == target_label)[0][0]
    # Step 2: train classifier
    clf = RandomForestClassifier(n_estimators=500, class_weight="balanced", n_jobs=-1, random_state=42)
    clf.fit(X, y_enc)
    # Step 3: get global feature importance
    importances = pd.Series(clf.feature_importances_, index=var_names)
    # Step 4: select top_k genes
    top_genes = importances.sort_values(ascending=False).head(top_k).index.tolist()
    # Step 5: optional filtering based on expression in other species same cell type
    species_specific_genes = []
    for gene in top_genes:
        expr_target = adatas[target_species][adatas[target_species].obs['Group']==target_cell, gene].layers["mean"]
        expr_other = []
        for sp, ad in adatas.items():
            if sp == target_species:
                continue
            mask = ad.obs['Group'] == target_cell
            if mask.sum() > 0:
                expr_other.append(ad[mask, gene].layers["mean"])
        if all(e <= expr_thresh for e in expr_other):
            species_specific_genes.append(gene)
    results = pd.DataFrame({
        "gene": top_genes,
        "importance": importances[top_genes].values,
        "is_species_specific": [g in species_specific_genes for g in top_genes]
    })
    ##
    return results

## 
wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs/"

##
anno_term = "Group"
expr_files = glob.glob(os.path.join(wkdir, anno_term.lower(), f"*_{anno_term.lower()}_mean_expr.h5ad"))

## Load expression data from hdf
adatas = {}
for file in expr_files:
    ##
    adata = ad.read_h5ad(file)
    ##
    species_name = re.search(r'/([^/]+)_'+anno_term.lower()+'_mean_expr\.h5ad$', file).group(1)
    adatas[species_name] = adata

##
hits = species_celltype_markers_filtered(adatas, target_species="human", target_cell="SN SOX6 Dopa")
