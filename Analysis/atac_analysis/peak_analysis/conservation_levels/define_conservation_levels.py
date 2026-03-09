import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

##
##
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Liftover results
liftover = pd.read_csv(os.path.join(analysis_dir, "annotations", "human_ref_liftover_HALPER_minMatch_0-5.tsv"), sep="\t")
liftover["region-id"] = liftover["human_ID"] + "-human"

## Read accessologs
accessologs = pd.read_csv(os.path.join(analysis_dir, "atac/accessologs", "accessologs_primate_orthologs.csv"))
accessologs.drop_duplicates(subset=['region'], inplace=True)

## Read in gini scores
gini_scores = pd.read_csv(os.path.join(analysis_dir, "atac/specificity", "gini_scores_combined.csv"))
gini_scores["region-id"] = (
    gini_scores["region"].astype(str) + "-" + gini_scores["species"].astype(str)
)

## Read in anndata with annotations
adata = ad.read_h5ad(os.path.join(analysis_dir, "cactus", "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"))

## ---------------------------
## Sequence conserved
## ---------------------------
ortholog_peaks = liftover.loc[liftover.ortholog == True,:].reset_index(drop=True)
ortholog_peaks.to_csv(os.path.join(analysis_dir, "atac/conservation", "sequence_conserved_peaks.csv"), index=False)

## ---------------------------
## Epi conserved
## ---------------------------
epi_conserved = accessologs.loc[accessologs.region.isin(ortholog_peaks["human_ID"]),:]
epi_conserved = epi_conserved.loc[epi_conserved.auroc > 0.7,:].reset_index(drop=True)
epi_conserved["human_ID"] = epi_conserved["region"].copy()

## Load in other species peak coordinates
for species in ["macaque", "marmoset"]:
    epi_conserved[f"{species}_ID"] = epi_conserved["human_ID"].map(
        dict(zip(
            liftover.loc[liftover[f"{species}_ID"].notna(), "human_ID"],
            liftover.loc[liftover[f"{species}_ID"].notna(), f"{species}_ID"]
        ))
    )
    epi_conserved[f"{species}_aligned_ID"] = epi_conserved["human_ID"].map(
        dict(zip(
            liftover.loc[liftover[f"{species}_ID"].notna(), "human_ID"],
            liftover.loc[liftover[f"{species}_ID"].notna(), f"{species}_aligned_ID"]
        ))
    )

## Load in gini scores
for species in ["human", "macaque", "marmoset"]:
    print(species)
    epi_conserved[f"{species}_gini"] = epi_conserved[f"{species}_ID"].map(
        dict(zip(
            gini_scores.loc[gini_scores.species == species, "region"],
            gini_scores.loc[gini_scores.species == species, "gini_scores"]
        ))
    )
epi_conserved["min_gini"] = epi_conserved[["human_gini", "macaque_gini", "marmoset_gini"]].min(axis=1)
epi_conserved["mean_gini"] = epi_conserved[["human_gini", "macaque_gini", "marmoset_gini"]].mean(axis=1)
epi_conserved.to_csv(os.path.join(analysis_dir, "atac/conservation", "epi_conserved_peaks.csv"), index=False)

## ---------------------------
## Epi conserved markers
## ---------------------------
epi_conserved_markers = epi_conserved.loc[epi_conserved["min_gini"] > 0.7,:].reset_index(drop=True)
epi_conserved_markers.to_csv(os.path.join(analysis_dir, "atac/conservation", "epi_conserved_marker_peaks.csv"), index=False)

## ---------------------------
## Species biased
## ---------------------------
species_biased_collection = {}
for species in ["human", "macaque", "marmoset"]:
    species_biased = epi_conserved.copy()
    species_biased["delta_gini"] = species_biased[f"{species}_gini"] - species_biased["min_gini"]
    species_biased = species_biased.sort_values(by=["delta_gini", f"{species}_gini"], ascending=False)
    species_biased = species_biased[
        (species_biased[f"{species}_gini"] > 0.6) & 
        (species_biased["delta_gini"] > 0.25)
    ].reset_index(drop=True)
    species_biased["species"] = species
    species_biased_collection[species] = species_biased

species_biased = pd.concat(species_biased_collection.values(), axis=0).reset_index(drop=True)
species_biased.to_csv(os.path.join(analysis_dir, "atac/conservation", "species_biased_peaks.csv"), index=False)

## ---------------------------
## Species specific
## ---------------------------
human_specific = liftover.loc[liftover.human_specific == True,:].reset_index(drop=True)
## Remove chrX,Y
human_specific = human_specific[~human_specific["Chromosome"].isin(["chrX", "chrY"])].reset_index(drop=True)
human_specific = human_specific[~human_specific["Chromosome"].str.contains("GL|KI")].reset_index(drop=True)

human_specific.to_csv(os.path.join(analysis_dir, "atac/conservation", "human_species_specific_peaks.csv"), index=False)