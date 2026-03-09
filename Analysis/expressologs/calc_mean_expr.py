import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import glob
import os

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/"

##### -- Load Data -- ######
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "Human_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "Marmoset_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
}

species_adatas = {}
for species in species_h5ads.keys(): 
    print(f"Loading {species} data...")
    adata = ad.read_h5ad(species_h5ads[species])
    adata.var_names_make_unique()
    species_adatas[species] = adata

## -------------------------------------------
## Loop through files computing mean expression
for anno_term in ["Group"]:
    for species in species_h5ads:
        ## Load species data
        print(f"Loading {species}")
        if os.path.exists(os.path.join(analysis_dir, anno_term.lower())) is False:
            os.mkdir(os.path.join(analysis_dir, anno_term.lower()))
        ## Compute mean expression per Subclass
        try:
            mean_expr = sc.get.aggregate(species_adatas[species], by=anno_term, func=["mean"])
            mean_expr.write_h5ad(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/mean_expr/{anno_term.lower()}/{species.lower()}_{anno_term.lower()}_mean_expr_orthologs.h5ad")
        except:
            print(f"Error in {species}")
