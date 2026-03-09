import os
from tqdm import tqdm
import pandas as pd
import anndata as ad
import anndata as ad
from memelite import fimo

## For plotting
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

species = "human"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/"
meme = os.path.join(analysis_dir, "annotations", "JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt")
fasta_file = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/{species}/peaks/merged_peaks.fasta"

## FIMO hits
hits = fimo(meme, fasta_file)
all_hits = pd.concat(hits)
all_hits["motif_TF"] = all_hits.motif_name.str.split(" ").str[1]

## Save hits (list of dataframes) as pkl 
all_hits.to_pickle(os.path.join(analysis_dir, "annotations", "TF_motifs", "fimo_hits.pkl"))

## Create motif by peak matrix
motif_mat = pd.DataFrame(0, index=all_hits.sequence_name.unique(), columns=all_hits.motif_TF.unique())
for _, row in tqdm(all_hits.iterrows()):
    motif_mat.at[row.sequence_name, row.motif_TF] += 1

##
adata = ad.AnnData(motif_mat.values)
adata.obs.index = motif_mat.index
adata.var.index = motif_mat.columns

## ---------------------------------------------
## Alignment stats and more!
## ---------------------------------------------
## Read in anndata with annotations
aligned_adata = ad.read_h5ad(os.path.join(analysis_dir, "cactus", "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"))

adata.obs = aligned_adata.obs.loc[adata.obs.index].copy()

## ---------------------------------------------
## TF class annotation via Lambert
## ---------------------------------------------
anno = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/TF_motifs/DatabaseExtract_v_1.01.xlsx", index_col=0)

## ---------------------------------------------
## Plot mean evo-distance and gini_scores per motif
## ---------------------------------------------
motif_stats = []
for motif in tqdm(adata.var.index):
    peaks_with_motif = adata.obs.index[adata[:, motif].X.flatten() > 0]
    mean_evo_distance = aligned_adata.obs.loc[peaks_with_motif, "evo_distance"].mean()
    mean_gini_score = aligned_adata.obs.loc[peaks_with_motif, "gini_scores"].mean()
    motif_stats.append({
        "motif_TF": motif,
        "mean_evo_distance": mean_evo_distance,
        "mean_gini_score": mean_gini_score,
        "num_peaks": len(peaks_with_motif)
    })

## Create a DataFrame from the motif_stats list
motif_stats_df = pd.DataFrame(motif_stats)
motif_stats_df["num_peaks_norm"] = motif_stats_df["num_peaks"] / aligned_adata.obs.shape[0]

##
motif_stats_df["anno"] = motif_stats_df["motif_TF"].str.upper().map(dict(zip(anno["HGNC symbol"], anno["DBD"])))
motif_stats_df["anno"] = motif_stats_df["anno"].fillna("Unknown")

## Set anno colors before plotting
anno_colors = {
    anno: to_hex(plt.cm.tab20(hash(anno) % 20)) if anno != "Unknown" else "#808080"
    for anno in motif_stats_df["anno"].unique()
}

## Set them into the dataframe after converting to RGB hex
motif_stats_df["anno_color"] = motif_stats_df["anno"].map(anno_colors)

motif_stats_df.to_csv(os.path.join(analysis_dir, "annotations", "TF_motifs", "motif_stats_with_annotations.csv"), index=False)

## ---------------------------------------------
## Plot
## ---------------------------------------------

## Plot the annotations in order of largest to smallest so the largest TF anno groups dont cover the smaller ones
anno_sizes = motif_stats_df.groupby("anno").size().sort_values(ascending=False).index.tolist()

plt.figure(figsize=(8, 8))
for anno in anno_sizes:
    subset = motif_stats_df[motif_stats_df["anno"] == anno]
    sc = plt.scatter(
        subset["mean_evo_distance"],
        subset["mean_gini_score"],
        s=subset["num_peaks_norm"] * 750,  # Scale size for visibility
        alpha=0.9 if anno != "Unknown" else 0.3,
        color=anno_colors[anno],
        label=anno
    )

plt.xlabel("Mean Evolutionary Distance")
plt.ylabel("Mean Gini Score")
plt.title("Motif Evolutionary Distance vs Gini Score")
plt.tight_layout()
plt.savefig(
    os.path.join(analysis_dir, "annotations", "TF_motifs", "motif_evo_distance_vs_gini_anno.pdf"),
    dpi=900
)



## Save the AnnData object
# adata.write_h5ad(f"{analysis_dir}/fimo_motif_by_peak.h5ad")
















# Add labels for TFs with high gini and low evolutionary distance
# highlight = motif_stats_df[
#     (motif_stats_df["mean_gini_score"] > 0.3) &
#     (motif_stats_df["mean_evo_distance"] < 0.025)
# ]

# for _, row in highlight.iterrows():
#     plt.annotate(
#         row["motif_TF"],  # adjust if your column is named differently
#         (row["mean_evo_distance"], row["mean_gini_score"]),
#         xytext=(5, 5),  # offset in points
#         textcoords="offset points",
#         fontsize=8,
#         color="black",
#         weight="bold"
#     )

# Add labels for TFs with high gini and low evolutionary distance
# highlight_low = motif_stats_df[
#     (motif_stats_df["mean_gini_score"] < 0.4) &
#     (motif_stats_df["mean_evo_distance"] > 0.035)
# ]

# highlight_low = motif_stats_df.loc[motif_stats_df["motif_TF"].str.contains("POU|SOX|HOX|EBF|SPI|OLIG|CTCF"), :]

# for _, row in highlight_low.iterrows():
#     plt.annotate(
#         row["motif_TF"],  # adjust if your column is named differently
#         (row["mean_evo_distance"], row["mean_gini_score"]),
#         xytext=(5, 5),  # offset in points
#         textcoords="offset points",
#         fontsize=2,
#         color="black",
#         weight="bold"
#     )