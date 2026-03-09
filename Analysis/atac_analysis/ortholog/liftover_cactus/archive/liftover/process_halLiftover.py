import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import gzip
import glob 
import os

import pyranges as pr
import pandas as pd

def summarize_halLiftover(halLiftover_df, species_peak_df, max_len=1000):
    """
    Use pyranges to summarize halLiftover results, keeping only the
    longest overlapping interval per ID.
    """
    ## Convert to PyRanges
    pr_all = pr.PyRanges(halLiftover_df)
    ## 
    results = []
    ## Process each ID separately
    for pid, sub in tqdm(pr_all.df.groupby("ID")):
        pr_sub = pr.PyRanges(sub)
        ## cluster overlapping intervals per chromosome
        clustered = pr_sub.cluster(slack=10)
        ## merge intervals within each cluster
        c_df = clustered.df
        merged = (
            c_df.groupby("Cluster", as_index=False, observed=True)
                .agg(
                    ID=("ID", "first"),
                    Chromosome=("Chromosome", "first"),
                    Start=("Start", "min"),
                    End=("End", "max")
                )
        )
        if merged.empty:
            continue
        ## pick the longest merged interval per chromosome while ignoring anything longer than 1kb
        merged["Length"] = merged.End - merged.Start
        longest = (
            merged[merged["Length"] <= max_len].sort_values("Length", ascending=False)
            .groupby("Chromosome", observed=True)
            .head(1)
        )
        longest = (
            merged[merged["Length"] <= max_len].sort_values("Length", ascending=False)
            .groupby("Chromosome", observed=True)
            .head(1)
        )
        ## keep the global longest across chromosomes for this ID
        best = longest.sort_values("Length", ascending=False).head(1)
        chrom, start, end, length = (
            best["Chromosome"].iloc[0],
            best["Start"].iloc[0],
            best["End"].iloc[0],
            best["Length"].iloc[0],
        )
        results.append((pid, chrom, length, start, end))
    ##
    per_chr = pd.DataFrame(results, columns=["ID", "best_chrom", "total_aligned_bp", "lifted_start", "lifted_end"])
    ## Merge with species peaks
    species_peak_df = species_peak_df.copy()
    species_peak_df["peak_length"] = species_peak_df["End"] - species_peak_df["Start"]
    ##
    result = species_peak_df.merge(per_chr, on="ID", how="left")
    result["total_aligned_bp"] = result["total_aligned_bp"].fillna(0).astype(int)
    result["best_chrom"] = result["best_chrom"].fillna("")
    result["lifted_start"] = result["lifted_start"].fillna(0).astype(int)
    result["lifted_end"] = result["lifted_end"].fillna(0).astype(int)
    result["frac_aligned"] = (
        result["total_aligned_bp"] / result["peak_length"].replace(0, pd.NA)
    ).fillna(0)
    result["total_aligned_bp_capped"] = result["total_aligned_bp"].clip(
        upper=result["peak_length"]
    )
    ##
    return result


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
species_liftovers = {
    "Mus_musculus": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/liftOver/Mus_musculusToHomo_sapiens.bed", 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
    "Macaca_mulatta": pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/liftOver/Macaca_mulattaToHomo_sapiens.bed", 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None)
}

##
# summarized = {}
# for species in species_liftovers.keys():
#     print(f"Processing {species}...")
#     res = summarize_halLiftover(species_liftovers[species], species_peaks[species])
#     ## We have ID to store the original species coordinates
#     res["Chromosome"] = res["best_chrom"]; del res["best_chrom"]
#     res["Start"] = res["lifted_start"]; del res["lifted_start"]
#     res["End"] = res["lifted_end"]; del res["lifted_end"]
#     summarized[species] = res
#     ##
#     res.to_csv(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/liftOver/{species}ToHomo_sapiens_processed.tsv", 
#                sep="\t", index=False)
    
## Read in results
summarized = {}
for species in species_liftovers.keys():
    print(f"Processing {species}...")
    summarized[species] = pd.read_csv(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/liftOver/{species}ToHomo_sapiens_processed.tsv", sep="\t")

## Now merge into liftover set anchoring on human when possible
def build_liftover_set(human_bed, species_liftovers, minMatch=0.5):
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
    human_df.set_index("ID", inplace=True, drop=False)
    ## Convert human peaks to PyRanges
    human_pr = pr.PyRanges(human_df)
    for species, df in species_liftovers.items():
        print(f"Processing {species}...")
        ## Remove rows without chromosome and with small overlap
        df = df[df["Chromosome"] != ""].copy()
        df = df[df.frac_aligned >= minMatch].copy() ## Only attempt merge in if >=50% of bases aligned to Human
        ## Assign alignment ID
        df["alignment_id"] = df["Chromosome"] + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
        ## Convert to pyrange
        species_pr = pr.PyRanges(df.rename(columns={"ID":f"{species}_ID"}))
        ## Find overlaps with human peaks
        intersect_df = human_pr.join(species_pr).as_df()
        human_df.loc[intersect_df.ID.tolist(), f"{species}_ID"] = intersect_df[f"{species}_ID"].values
        human_df.loc[intersect_df.ID.tolist(), f"{species}_aligned_ID"] = intersect_df["alignment_id"].values
    ## Rename human ID column
    human_df.rename(columns={"ID":"human_ID"}, inplace=True)
    ## Return
    return human_df

## Create human-anchored liftover set
liftover_set = build_liftover_set(
    species_peaks["Homo_sapiens"],
    summarized,
    minMatch=0.5
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

## Save for CACTUS analysis
liftover_set.to_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/annotations/human_ref_liftover_minMatch_0-5.tsv", sep="\t", index=False, header=True)
