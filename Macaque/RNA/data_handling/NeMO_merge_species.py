import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import os

## ----------------------------------------
cross_species = {"Human": ad.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/Human_basalganglia_NeMO_BICAN.h5ad"),
                 "Macaque": ad.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia/Macaque_basalganglia_NeMO_BICAN.h5ad")}


## Loop through species and update gene names to orthologs
ortholog_table = pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes/all_gene_ids.csv', header=0)
for species, species_name in zip(['Macaque', 'Human'], 
                                 ['Macaca_nemestrina', "Homo_sapiens"]):
    ##
    print("Identifying orthologs for: " + species)
    ortholog_dict = dict(zip(ortholog_table.loc[ortholog_table['species']==species_name,'gene'], 
                             ortholog_table.loc[ortholog_table['species']==species_name,'ortholog_symbol']))
    cross_species[species].var['mgi_symbol'] = cross_species[species].var.index
    cross_species[species].var.index = cross_species[species].var.index.to_series().replace(ortholog_dict)
    cross_species[species].var_names_make_unique()

##
adatas = list(cross_species.values())
adata_merged = ad.concat(adatas, merge="same", uns_merge="same") 

##
adata_merged.write_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/basalganglia_xspecies_HMBA.h5ad")