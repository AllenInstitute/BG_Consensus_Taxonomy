import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/ATAC-conservation/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"

## spinalcord metadata
cluster_meta = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Identify species-biased peaks and plot heatmaps
top_biased_hits = {}
for ref_species in ["human", "macaque", "mouse"]:
    top_biased_hits[ref_species] = pd.read_csv(os.path.join(analysis_dir, "specificity", f"{ref_species}_biased_peaks.csv"), index_col=0)

all_top_biased = pd.concat(top_biased_hits, axis=0).reset_index(drop=True)

## Level 3 hits
level3_hits = pd.read_csv(os.path.join(analysis_dir, "specificity", "level3_ortholog_peaks_per_cell_type_human.csv"), index_col=0)

## I want to see the counts of peaks per cell type per species in the biased peaks
level3_counts = level3_hits.groupby("cell_type").size().reset_index(name="num_peaks")

## For the species biased peaks
biased_counts = all_top_biased.groupby(["species", "max_accessible_celltype"]).size().reset_index(name="num_peaks")

## Now combined into one dataframe with cell types as rows and species/level3 as columns
combined_counts = biased_counts.pivot(index="max_accessible_celltype", columns="species", values="num_peaks").fillna(0).reset_index()
combined_counts = combined_counts.merge(level3_counts, left_on="max_accessible_celltype", right_on="cell_type", how="left").rename(columns={"num_peaks": "level3_num_peaks"}).drop(columns=["cell_type"]).fillna(0) 

## Save combined counts
combined_counts.to_csv(os.path.join(data_dir, "peak_counts_by_celltype_and_species_biased_level3.csv"))