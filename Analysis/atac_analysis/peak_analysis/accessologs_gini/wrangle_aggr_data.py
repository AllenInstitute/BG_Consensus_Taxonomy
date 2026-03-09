import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import glob
import os

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data"

## -------------------------------------------
## Identify species files, could also be a list from glob
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "crested_adata", "human_spinalcord_hmba_crested.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "crested_adata", "macaque_spinalcord_hmba_crested.h5ad"),
    "mouse": os.path.join(base_dir, "mouse", "crested_adata", "mouse_spinalcord_hmba_crested.h5ad"),
}

## -------------------------------------------
## Define the ortholog peak ids
peakDir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks"
peak_lookup = pd.read_csv(peakDir + f"/peak_lookup_ref_human.csv")

## Macaque chrom alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## -------------------------------------------
## Loop through files computing mean expression
for species in species_h5ads:
    ## Load species data
    print(f"Loading {species}")
    adata = ad.read_h5ad(species_h5ads[species])
    if species == "human":
        adata.var["peak_id"] = adata.var.index
        continue
    elif species == "macaque":
        adata.var["chr"].replace(dict(zip(macaque_chrom_alias[2], macaque_chrom_alias[0])), inplace=True)
        adata.var["peak_id"] = adata.var["chr"].astype(str) + ":" + adata.var["start"].astype(str) + "-" + adata.var["end"].astype(str)
        adata.var["peak_id"] = adata.var["peak_id"].map(peak_lookup.set_index("human_peak_name")[species])
    else:
        adata.var["peak_id"] = adata.var.index
        adata.var["peak_id"] = adata.var["peak_id"].map(peak_lookup.set_index("human_peak_name")[species])
    ## Drop NaN var_names
    adata.var.dropna(subset=["peak_id"], inplace=True)
    ##
    adata.write_h5ad(peakDir + f"/{species}_spinalcord_hmba_crested_accessologs.h5ad")

    # ## Compute mean expression per Subclass
    # try:
    #     mean_expr = sc.get.aggregate(adata, by=anno_term, func=["mean"])
    #     mean_expr.write_h5ad(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Expressologs/{species.lower()}_{anno_term.lower()}_mean_expr.h5ad")
    # except:
    #     print(f"Error in {species}")
