import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import umap
import matplotlib.pyplot as plt
import gzip
import glob 
import os
import anndata as ad
from scipy.sparse import csr_matrix
from ete3 import NCBITaxa

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"

## 
species_peak_file = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed"
path_to_species_overlaps = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/archive/analysis"

## Read in overlap data
adata = ad.read_h5ad(os.path.join(path_to_species_overlaps, "human_cactus_alignment_basepair_overlap_capped.h5ad"))

## Pull taxa from NCBI
ncbi = NCBITaxa()

## Add in taxanomic info for each species
species_list = list(adata.var.index)  # e.g., 'Homo_sapiens', 'Mus_musculus'
taxid_dict = {}

adata.var.index = adata.var.index.str.replace("merged_peaks.", "")
adata.var.index = adata.var.index.str.replace(".HALPER.narrowPeak.gz", "")
adata.var.index = adata.var.index.str.replace("_", " ")

## Gather taxa info
for sp in tqdm(species_list):
    try:
        taxid = ncbi.get_name_translator([sp.replace("_", " ")])[sp.replace("_", " ")][0]
        lineage = ncbi.get_lineage(taxid)              # list of taxids up to root
        names = ncbi.get_taxid_translator(lineage)     # map taxid -> name
        ranks = ncbi.get_rank(lineage)                 # taxid -> rank
        taxid_dict[sp] = {ranks[t]: names[t] for t in lineage}
    except:
        taxid_dict[sp] = {}

## Convert to DataFrame
taxa_df = pd.DataFrame.from_dict(taxid_dict, orient='index')

## Remove no rank column
taxa_df = taxa_df.drop(columns=["no rank", "cellular root"])

## Add taxa info
adata.var = pd.concat([adata.var, taxa_df], axis=1)

adata.write_h5ad(os.path.join(path_to_species_overlaps, "human_cactus_alignment_basepair_overlap_capped.h5ad"))