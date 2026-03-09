import os
import sys
import crested
import keras
import glob
import re
import pandas as pd
import pyranges as pr
import pysam
import anndata as ad
import tensorflow as tf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse
tf.config.list_physical_devices('GPU')

## Function to pad sequences equally on both sides
def pad_sequence(sequence, target_width):
    ## If sequence is to long cut out the middle
    if len(sequence) >= target_width:
        center = len(sequence) // 2
        start = (centers - (target_width / 2)).astype(int)
        end = (centers + (target_width / 2)).astype(int)
        return sequence[start:end]
    ## Otherwise pad the sequence with N's
    total_padding = target_width - len(sequence)
    left_padding = total_padding // 2
    right_padding = total_padding - left_padding
    return ('N' * left_padding) + sequence + ('N' * right_padding)

model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"

## -------------------------
enhancer_table = pd.read_excel("/home/nelson.johansen/Analysis/HMBA_Genomics/Analysis/SeqModel/enhancer_tables/striatal_tools_supp_table1.xlsx")
enhancer_table["Genome assembly"] = enhancer_table["Genome assembly"].str.lower()

## Clean up enhancer IDs
enhancer_table["enhancer ID"] = enhancer_table["enhancer ID"].str.strip()
enhancer_table["enhancer_type"] = "full length"
enhancer_table.loc[enhancer_table["enhancer ID"].str.contains("_"), "enhancer_type"] = "optimized"

## Sequence handling
enhancer_table["Sequence"] = enhancer_table["Enhancer sequence"]
enhancer_table = enhancer_table.loc[enhancer_table["Sequence"].notna()]

## --------------------------
top_peaks_per_type = pd.read_csv("/home/nelson.johansen/Analysis/HMBA_Genomics/Analysis/SeqModel/enhancer_tables/top_peaks_per_type.csv")

## Remove "," from location column
# enhancer_table["location"] = enhancer_table["location"].str.replace(",", "", regex=False)

## Extract chromosome, start, and end from the "Region" column
# enhancer_table["Chromosome"] = enhancer_table["location"].str.split(":").str[0]
# enhancer_table["Start"] = enhancer_table["location"].str.split(":").str[1].str.split("-").str[0].astype(int)
# enhancer_table["End"] = enhancer_table["location"].str.split(":").str[1].str.split("-").str[1].astype(int)

## Gather genomes
# genomes_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes"
# fasta_dict = {
#     "hg38": pysam.FastaFile(os.path.join(genomes_dir, "human", "fasta", "hg38.fa")),
#     "mm10": pysam.FastaFile(os.path.join(genomes_dir, "mouse", "fasta", "mm10.fa")),
#     "rhemac10": pysam.FastaFile(os.path.join(genomes_dir, "macaque", "fasta", "genome.fa")),
# }

## Pull sequence
# enhancer_table["Sequence"] = enhancer_table.apply(
#     lambda row: fasta_dict[row["Genome build"]].fetch(row["Chromosome"], row["Start"], row["End"]),
#     axis=1
# )

## -------------------------
# Convert all sequences to CREsted model length
enhancer_table["Sequence_CRESted"] = enhancer_table["Sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
model = tf.keras.models.load_model(os.path.join(model_dir, "macaque/finetune/11.keras"), compile=False)

## -------------------------
model_embedding = crested.tl.extract_layer_embeddings(
    input = enhancer_table["Sequence_CRESted"].tolist(),
    model = model,
    layer_name = "global_average_pooling1d",
)

## --------------------------
# Plot UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

# --- Parameters ---
method = "umap" 

# --- Dimensionality reduction ---
if method == "umap":
    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(model_embedding)
elif method == "tsne":
    reducer = TSNE(n_components=2, random_state=42, perplexity=30)
    embedding = reducer.fit_transform(model_embedding)
else:
    raise ValueError("Invalid method. Choose 'umap' or 'tsne'.")

## --- Build DataFrame ---
embedding_df = pd.DataFrame(embedding, columns=[f"{method.upper()}1", f"{method.upper()}2"])
embedding_df["target_celltype"] = enhancer_table["Cell type"]
embedding_df["enhancer_type"] = enhancer_table["enhancer_type"]

## convert nan to "Unknown""
embedding_df.to_csv(f"/home/nelson.johansen/Striatal_tools_Enhancer_Set_{method}_embedding.csv", index=False)

## Create a color palette (e.g., with seaborn)
labels = embedding_df["target_celltype"].unique()
palette = sns.color_palette("tab20", n_colors=len(labels))
color_map = dict(zip(labels, palette))

# --- Plotting ---
plt.figure(figsize=(9, 6), dpi=300)
for celltype in list(color_map.keys())[::-1]:
    print(f"Plotting cell type: {celltype}")
    subset = embedding_df[embedding_df["target_celltype"] == celltype]
    plt.scatter(
        subset[f"{method.upper()}1"],
        subset[f"{method.upper()}2"],
        color=color_map[celltype],
        alpha=0.6,
        rasterized=True,
        s=20,
        label=celltype,
    )

plt.title(f"{method.upper()} Visualization of CRESted Model Embeddings")
plt.xlabel(f"{method.upper()}1")
plt.ylabel(f"{method.upper()}2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Target population")
plt.tight_layout()
plt.savefig('/home/nelson.johansen/seq_umap_striatal_tools.pdf', bbox_inches='tight')
