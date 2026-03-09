
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import glob 
from functools import reduce
import re
import os
from collections import Counter
import scvi

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/samap"

##
cross_species_paths = {
    "HMBA:Human" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/human/Human_HMBA_basalganglia_AIT_pre-print.h5ad",
    "HMBA:Macaque" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/Macaque_HMBA_basalganglia_AIT_pre-print.h5ad",
    "HMBA:Marmoset": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/marmoset/Marmoset_HMBA_basalganglia_AIT_pre-print.h5ad",
    "Mouse_Broad": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/mouse/Broad_Mouse_BasalGanglia_pre-print.h5ad"
}

## ----------------------------------------
## Load cross species data
cross_species = {}
for k in cross_species_paths.keys():
    adata = ad.read_h5ad(cross_species_paths[k])
    if k == "Mouse_Broad":
        adata.var_names = adata.var.gene_name
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        adata.layers["UMIs"] = adata.X.copy()  # Store the raw counts in a layer
        adata.obs["donor_id"] = "mouse"
        adata.obs["organism_sci"] = "Mus musculus"
        adata = adata[adata.obs["AIT21.subclass"].isin(["055 STR Lhx8 Gaba"])].copy()
        del adata.X
    else:
        adata.layers["UMIs"] = adata.raw.X.copy()  # Store the raw counts in a layer
        adata = adata[adata.obs.Group.isin(["STR-BF TAC3-PLPP4-LHX8 GABA", "STR TAC3-PLPP4 GABA"])].copy()
        del adata.raw
        del adata.X
    cross_species[k] = adata

for k in cross_species.keys():
    print(f"{k}: {cross_species[k].shape}")
    ## ----------------------------------------
    cross_species[k].write_h5ad(os.path.join(analysis_dir, f"{k.replace(':','_')}_for_samap.h5ad"))