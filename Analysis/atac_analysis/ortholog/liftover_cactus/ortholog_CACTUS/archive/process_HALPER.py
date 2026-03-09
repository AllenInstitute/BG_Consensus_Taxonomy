import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import gzip
import glob 
import os

def read_bed_gz(path):
    with gzip.open(path, "rt") as f:
        return pd.read_csv(f, sep="\t", header=None)
    
## Read in command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--species_id", type=str, help="Index of species to process", default="0")
parser.add_argument("--ref_species", type=str, help="Reference species, e.g. human", default="human")
args = parser.parse_args()

## Paths to species peak files
species_peak_file = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/{args.ref_species}/Consensus_peaks/merged_peaks_with_names.bed"

## Load reference species peaks
species_bed = pd.read_csv(species_peak_file, sep="\t", header=None)
species_bed.columns = ["Chromosome", "Start", "End", "ID"]
species_bed["peak_length"] = species_bed["End"] - species_bed["Start"] - 1

## Gather all tFiles for this species
path_to_species_tfile = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/{args.ref_species}"
speciesTo_tFiles = [f for f in glob.glob(f"{path_to_species_tfile}/*tFile.bed.gz") if "Anc" not in f]

## Make analysis directory if it doesn't exist
os.makedirs(os.path.join(path_to_species_tfile, "analysis"), exist_ok=True)

##
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
        ## If no alignment then continue
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
        ## If no best (e.g. all > max_len), set to 0
        if best.empty:
            continue
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

## for speciesTo_tFile in tqdm(speciesTo_tFiles):
speciesTo_tFile = speciesTo_tFiles[int(args.species_id)]
species = speciesTo_tFile.split("/")[-1].split(".")[1]

##
print(f"Processing {species} for {args.ref_species}...")

## Read in halLiftover results
halLiftover_df = read_bed_gz(speciesTo_tFile)
halLiftover_df.columns = ["Chromosome", "Start", "End", "ID", "Score", "Strand"]

## Compute overlaps and summarize
result = summarize_halLiftover(halLiftover_df, species_bed)

## Save results
result.to_csv(os.path.join(path_to_species_tfile, "analysis", f"{species}_human_peaks_alignment_summary.tsv"), sep="\t", index=False)




##
# def summarize_halLiftover(halLiftover_df, species_peak_df):
#     """
#     Use pyranges to summarize halLiftover results, keeping only the
#     longest overlapping interval per ID.
#     """
#     ## Convert to PyRanges
#     pr_all = pr.PyRanges(halLiftover_df)
#     ## 
#     results = []
#     ## Process each ID separately
#     for pid, sub in tqdm(pr_all.df.groupby("ID")):
#         pr_sub = pr.PyRanges(sub)
#         ## cluster overlapping intervals per chromosome
#         clustered = pr_sub.cluster(slack=10)
#         ## merge intervals within each cluster
#         c_df = clustered.df
#         merged = (
#             c_df.groupby("Cluster", as_index=False, observed=True)
#                 .agg(
#                     ID=("ID", "first"),
#                     Chromosome=("Chromosome", "first"),
#                     Start=("Start", "min"),
#                     End=("End", "max")
#                 )
#         )
#         if merged.empty:
#             continue
#         ## pick the longest merged interval per chromosome
#         merged["Length"] = merged.End - merged.Start
#         longest = (
#             merged.sort_values("Length", ascending=False)
#             .groupby("Chromosome", observed=True)
#             .head(1)
#         )
#         longest = (
#             merged.sort_values("Length", ascending=False)
#             .groupby("Chromosome", observed=True)
#             .head(1)
#         )
#         ## keep the global longest across chromosomes for this ID
#         best = longest.sort_values("Length", ascending=False).head(1)
#         chrom, start, end, length = (
#             best["Chromosome"].iloc[0],
#             best["Start"].iloc[0],
#             best["End"].iloc[0],
#             best["Length"].iloc[0],
#         )
#         results.append((pid, chrom, length, start, end))
#     ##
#     per_chr = pd.DataFrame(results, columns=["ID", "best_chrom", "total_aligned_bp", "lifted_start", "lifted_end"])
#     ## Merge with species peaks
#     species_peak_df = species_peak_df.copy()
#     species_peak_df["peak_length"] = species_peak_df["End"] - species_peak_df["Start"]
#     ##
#     result = species_peak_df.merge(per_chr, on="ID", how="left")
#     result["total_aligned_bp"] = result["total_aligned_bp"].fillna(0).astype(int)
#     result["best_chrom"] = result["best_chrom"].fillna("")
#     result["lifted_start"] = result["lifted_start"].fillna(0).astype(int)
#     result["lifted_end"] = result["lifted_end"].fillna(0).astype(int)
#     result["frac_aligned"] = (
#         result["total_aligned_bp"] / result["peak_length"].replace(0, pd.NA)
#     ).fillna(0)
#     result["total_aligned_bp_capped"] = result["total_aligned_bp"].clip(
#         upper=result["peak_length"]
#     )
#     ##
#     return result