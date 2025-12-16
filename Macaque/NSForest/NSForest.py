import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.stats import entropy
import os, sys

import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating as ev
from nsforest import plotting as pl
import os

##
def grep(l, s):
    return [i for i in l if s in i]

##
species = "Macaque"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

## Save donor corrected latent space
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_HMBA_AIT11-9_anno_031525.h5ad"))

## Compute "conserved" marker genes for each anno level
for anno in ["cluster_id"]:
    ## Compute medians of `anno` per species
    adata = pp.prep_medians(adata, anno, use_mean = False, positive_genes_only = True)
    ## Compute beta scores grouping on "species"
    adata = pp.prep_binary_scores(adata, anno, medians_header = "medians_")
    ## Create folder for outputs
    output_path = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NSForest/{anno}/"
    os.makedirs(output_path, exist_ok=True)
    ##
    results = nsforesting.NSForest(adata, 
                                    cluster_header = anno, 
                                    save_supplementary = True, 
                                    output_folder = output_path, 
                                    outputfilename_prefix = f"{species}_BG_{anno}_NSForest")



# ## Compute entropy on beta scores to identify Class-specific markers across species
# def compute_row_entropy(row):
#     # Normalize the row values to probabilities
#     probabilities = row / row.sum()
#     # Use scipy.stats.entropy (base 2 for Shannon entropy)
#     return entropy(probabilities, base=2)
# # Compute entropy for each row
# gene_binary_score_entropies = adata.varm["binary_scores_Subclass"].apply(compute_row_entropy, axis=1)


# marker_genes_dict = {
#     "Cholinergic": ["CHAT", "NOS1"],
#     "Motorneurons": ["CHAT", "PRPH"],
#     "GABAergic": ["GAD1", "SLC6A5"],
#     "Glutamatergic": ["SLC17A6", "CPNE4"],
#     "Non-Neurons": ["ZEB2"],
# }

# from matplotlib import rc_context
# with rc_context({"figure.dpi": 600}):
#     sc.pl.dotplot(adata, marker_genes_dict, "Class", dendrogram=False, save="Class_markers.png")             