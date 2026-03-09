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
parser.add_argument("--species_id", type=str, help="Index of species to process")

## 
species_peak_file = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed"
path_to_species_tfile = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/archive"
speciesTo_tFiles = [f for f in glob.glob(f"{path_to_species_tfile}/*tFile.bed.gz") if "Anc" not in f]

# Load original human peaks
human_bed = pd.read_csv(species_peak_file, sep="\t", header=None)
human_bed.columns = ["chrom", "start", "end", "id"]
human_bed["peak_length"] = human_bed["end"] - human_bed["start"] - 1

##
def summarize_halLiftover(halLiftover_df, species_peak_df):
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
        ## pick the longest merged interval per chromosome
        merged["Length"] = merged.End - merged.Start
        longest = (
            merged.sort_values("Length", ascending=False)
            .groupby("Chromosome", observed=True)
            .head(1)
        )
        longest = (
            merged.sort_values("Length", ascending=False)
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

## for speciesTo_tFile in tqdm(speciesTo_tFiles):
speciesTo_tFile = speciesTo_tFiles[int(parser.parse_args().species_id)]
species = speciesTo_tFile.split("/")[-1].split(".")[1]

##
print(f"Processing {species} ...")
result = summarize_halLiftover(speciesTo_tFile, human_bed)
result.to_csv(os.path.join(path_to_species_tfile, "analysis", f"{species}_human_peaks_alignment_summary.tsv"), sep="\t", index=False)

## Combine all species results from tFile_dfs by creating a matrix with peaks x species
# results = pd.DataFrame()
# results.index = human_bed["id"]
# for species, df in tqdm(tFile_dfs.items()):
#     results[f"{species}"] = df["total_aligned_bp_capped"].values

# ## Create a UMAP
# import umap
# reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean')
# embedding = reducer.fit_transform(results.fillna(0).values)
# embedding_df = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"], index=results.index)

# ## Plot UMAP
# import matplotlib.pyplot as plt
# plt.figure(figsize=(8,6))
# plt.scatter(embedding_df["UMAP1"], embedding_df["UMAP2"], s=1, alpha=0.5)
# plt.title("UMAP of Human Peaks based on Alignment to Other Species")
# plt.xlabel("UMAP1")
# plt.ylabel("UMAP2")
# plt.savefig("/home/nelson.johansen/human_peaks_umap.png", dpi=300)

# ## Calculate the number of elements per species with >= 90% alignment and <= 10% alignment
# summary_stats = {}
# for species, df in tFile_dfs.items():
#     n_peaks_90 = (df["frac_aligned"] >= 0.9).sum()
#     n_peaks_10 = (df["frac_aligned"] <= 0.1).sum()
#     summary_stats[species] = {
#         "n_peaks_>=90%_aligned": n_peaks_90,
#         "n_peaks_<=10%_aligned": n_peaks_10,
#         "total_peaks": len(df)
#     }

# # Save output
# result.to_csv("human_peaks_alignment_summary.tsv", sep="\t", index=False)
