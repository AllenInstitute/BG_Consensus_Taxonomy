import sys, os
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py

## Helpful locations which are assumed to already exist
annodir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
datadir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"

##
adata_obs = pd.read_csv(os.path.join(datadir, "xspecies", "consensus_hmba_basalganglia_metadata.csv"), index_col=0)

## Abundance 
cell_counts = adata.obs.groupby(['Group', 'load_id']).size().reset_index(name='n_cells_group_per_library')

## Add species based on donor_id
species_meta = adata.obs.drop_duplicates(subset=['load_id', 'Group']).copy()
cell_counts['organism'] = cell_counts['load_id'].replace(dict(zip(
    species_meta.load_id,
   species_meta.organism)))
cell_counts['roi'] = cell_counts['load_id'].replace(dict(zip(
    species_meta.load_id,
   species_meta.anatomical_region)))
cell_counts = cell_counts[["Group", "load_id", "organism", "n_cells_group_per_library"]]

## Grep for "Dopa" in "Group" to filter down to Dopa groups
cell_counts = cell_counts[cell_counts['Group'].str.contains("Dopa")].copy()

## Add species counts
cell_counts_species = adata.obs.groupby(['Group', 'organism']).size().reset_index(name='n_cells_group_per_organism')
cell_counts_species = adata.obs.groupby(['load_id']).size().reset_index(name='n_cells_library')

## Save
cell_counts.to_csv("/home/nelson.johansen/cell_counts_per_library.csv", index=False)
cell_counts_species.to_csv("/home/nelson.johansen/cell_counts_per_species.csv", index=False)
