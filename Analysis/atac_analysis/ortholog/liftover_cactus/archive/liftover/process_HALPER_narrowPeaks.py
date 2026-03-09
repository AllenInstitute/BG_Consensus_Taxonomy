import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import gzip
import glob 
import os

import pyranges as pr
import pandas as pd

## Peak data
species_peaks = {
    "Mus_musculus": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/mouse/Consensus_peaks/merged_peaks_with_names.bed", 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
    "Macaca_mulatta": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/macaque/Consensus_peaks/merged_peaks_with_names.bed", 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
    "Homo_sapiens": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed", 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
}

## halLiftover data
names = ["Chromosome", "Start", "End", "summit_pos", "ID", "target_length", "query_length", "target_summit_to_target_ortholog_start_length", "target_summit_to_target_ortholog_end_length"]
species_liftovers = {
    "Mus_musculus": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/mouse/merged_peaks_with_names.Mus_musculusToHomo_sapiens.HALPER_2.narrowPeak.gz", 
                                    names=names, sep="\t", header=None, compression="gzip"),
    "Macaca_mulatta": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/macaque/merged_peaks_with_names.Macaca_mulattaToHomo_sapiens.HALPER_2.narrowPeak.gz", 
                                  names=names, sep="\t", header=None, compression="gzip")
}

## Now merge into liftover set anchoring on human when possible
def build_liftover_set(human_bed, species_liftovers, minMatch=250):
    """
    Build human-anchored liftover peak set.

    Parameters
    ----------
    human_bed : pd.DataFrame
        Human peaks with columns ["Chromosome","Start","End","ID"].
    species_liftovers : dict[str, pd.DataFrame]
        Dict mapping species -> liftover peaks in human coords
        (each with ["Chromosome","Start","End","ID"]).

    Returns
    -------
    pd.DataFrame
        Ortholog peak set with columns:
        ["ortholog_region_id","Chromosome","Start","End","source","ID"]
    """
    ## Anchor on human peaks
    human_df = human_bed.copy()
    ## Set index without giving a name
    human_df.set_index("ID", inplace=True, drop=False)
    ## Convert human peaks to PyRanges
    human_pr = pr.PyRanges(human_df)
    for species, df in species_liftovers.items():
        print(f"Processing {species}...")
        ## Remove rows without chromosome and with small overlap
        df = df[df["Chromosome"] != ""].copy()
        df = df[df.target_length >= minMatch].copy() ## Only attempt merge in if >=50% of bases aligned to Human
        ## Assign alignment ID
        df["alignment_id"] = df["Chromosome"] + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
        ## Convert to pyrange
        species_pr = pr.PyRanges(df.rename(columns={"ID":f"{species}_ID"}))
        ## Find overlaps with human peaks
        intersect_df = human_pr.join(species_pr).as_df()
        #intersect_df = human_pr.join(species_pr, how="right").as_df()
        #intersect_df = intersect_df[intersect_df["ID"] != "-1"].copy()
        human_df.loc[intersect_df.ID.tolist(), f"{species}_ID"] = intersect_df[f"{species}_ID"].values
        human_df.loc[intersect_df.ID.tolist(), f"{species}_aligned_ID"] = intersect_df["alignment_id"].values
    ## Rename human ID column
    human_df.rename(columns={"ID":"human_ID"}, inplace=True)
    ## Return
    return human_df

## Create human-anchored liftover set
liftover_set = build_liftover_set(
    species_peaks["Homo_sapiens"],
    species_liftovers,
    minMatch=250 ## 50% overlap for 500bp peaks
)

## Identify human specific peaks by is.na on all but human_ID columns
liftover_set["human_specific"] = False
liftover_set.loc[
    (liftover_set["Mus_musculus_ID"].isna()) & 
    (liftover_set["Macaca_mulatta_ID"].isna()),
    "human_specific"
] = True

liftover_set["ortholog"] = False
liftover_set.loc[
    (~liftover_set["Mus_musculus_ID"].isna()) & 
    (~liftover_set["Macaca_mulatta_ID"].isna()),
    "ortholog"
] = True

liftover_set.columns = liftover_set.columns.str.replace("Mus_musculus", "mouse", regex=False)
liftover_set.columns = liftover_set.columns.str.replace("Macaca_mulatta", "macaque", regex=False)

## ortholog_set
ortholog_set = liftover_set[liftover_set["ortholog"] == True]
len(ortholog_set)

## Save for CACTUS analysis
liftover_set.to_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/annotations/human_ref_liftover_HALPER_minMatch_0-5.tsv", sep="\t", index=False, header=True)