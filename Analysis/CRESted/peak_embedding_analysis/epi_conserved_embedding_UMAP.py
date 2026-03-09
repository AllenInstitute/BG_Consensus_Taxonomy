import os
import sys
import glob
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import sparse
import tensorflow as tf
import keras
import crested
import anndata as ad
import pyranges as pr
import pysam
from pyfaidx import Fasta

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

# === TensorFlow Setup ===
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

## Handle nans
def parse_coordinates(coord):
    if pd.isna(coord) or not isinstance(coord, str) or ":" not in coord or "-" not in coord:
        return pd.Series([np.nan, np.nan, np.nan])
    try:
        chrom, pos = coord.split(":")
        start, end = pos.split("-")
        return pd.Series([chrom, int(start), int(end)])
    except Exception:
        return pd.Series([np.nan, np.nan, np.nan])

# ## Function to fetch sequences in batch
# def fetch_sequences(chrom, starts, ends):
#     ref = fasta[chrom]
#     return [ref[start:end].seq for start, end in zip(starts, ends)]

# ## Function to fetch sequences in order
# def fetch_sequences(gr_df, fasta, fetch_fn, seq_len=500):
#     """
#     Fetch sequences from a FASTA file for regions in a PyRanges dataframe,
#     preserving the original order.

#     Parameters
#     ----------
#     gr_df : pandas.DataFrame
#         A dataframe with columns ['Chromosome', 'Start', 'End'] at minimum.
#     fasta : Fasta
#         An open pyfaidx.Fasta object.
#     fetch_fn : function
#         A function with signature (chrom, starts, ends) -> list of sequences.
#     seq_len : int, optional
#         Length of sequence to insert if chromosome not found. Default 500.

#     Returns
#     -------
#     list of str
#         Sequences in the same order as the original gr_df.
#     """
#     gr_df["orig_idx"] = gr_df.index
#     seqs_per_idx = {}
#     for chrom, subdf in tqdm(gr_df.groupby("Chromosome"), desc="Fetching sequences"):
#         if chrom not in fasta.keys():
#             seqs_chr = ["N" * seq_len] * len(subdf)
#         else:
#             seqs_chr = fetch_fn(chrom, subdf["Start"], subdf["End"])
#         seqs_per_idx.update(dict(zip(subdf["orig_idx"], seqs_chr)))
#     ## Reorder to match original index
#     ordered_seqs = [seqs_per_idx[i] for i in sorted(gr_df["orig_idx"])]
#     del gr_df["orig_idx"]
#     return ordered_seqs

## From CRESted
def calc_gini(targets: np.ndarray) -> np.ndarray:
    """
    Return Gini scores for the given targets.

    This function calculates the Gini coefficient for each row in the targets array
    and assigns the score to the maximum value's position in the corresponding row
    of the Gini scores array.

    Parameters
    ----------
    targets
        A 2D numpy array where each row represents a set of target values.

    Returns
    -------
    gini scores
        A 2D numpy array with the same shape as `targets` containing Gini scores,
        where each score is assigned to the position of the maximum value in each row.
    """
    def _gini(array: np.ndarray) -> float:
        """Calculate the Gini coefficient of a numpy array."""
        array = (
            array.flatten().clip(0, None) + 0.0000001
        )  # Ensure non-negative values and avoid zero
        array = np.sort(array)
        index = np.arange(1, array.size + 1)
        return (np.sum((2 * index - array.size - 1) * array)) / (
            array.size * np.sum(array)
        )
    ##
    gini_scores = np.zeros_like(targets)
    ##
    for region_idx in range(targets.shape[0]):
        region_scores = targets[region_idx]
        max_idx = np.argmax(region_scores)
        gini_scores[region_idx, max_idx] = _gini(region_scores)
    ##
    return gini_scores

## -------------------------
##
## -------------------------
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
## Load data
## -------------------------
adata = ad.read_h5ad(os.path.join(data_dir, "human/crested_adata/human_basalganglia_hmba_pre-print_crested.h5ad"))

## -------------------------
## Genome
## -------------------------
species = "human"
species_to_reference = {
    "human": "hg38",
    "mouse": "mm10",
    "macaque": "rheMac10",
    "marmoset": "calJac4"
}
## -------------------------
reference = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/{species}"
if species == "macaque":
    ## Macaque specific reference
    reference = os.path.join(reference, "ncbi")
    ##
    chrom_alias = pd.read_table(os.path.join(reference, "rheMac10.chromAlias.txt"), delimiter="\t")
    chrom_alias.columns = ["seqName", "aliasName", "UCSC"]
    ## Gather chr1-20, chrX,Y
    chromosomes = chrom_alias.loc[chrom_alias.UCSC.str.startswith("NC_")]
    chromosomes = chromosomes.loc[~chromosomes.seqName.str.contains("chrM")]
elif species == "marmoset":
    reference = os.path.join(reference, "hmba")

## Fasta and GTF paths
fasta_path = os.path.abspath(os.path.join(reference, f"fasta/genome.fa"))
gtf_path = os.path.abspath(os.path.join(reference, f"genes/genes.gtf.gz"))

## Load genome FASTA
fasta = Fasta(fasta_path)

genome = crested.Genome(
    fasta=fasta_path,
    annotation=gtf_path,
)

## -------------------------
## Load in enhancer table
## -------------------------
enhancer_table = pd.read_excel(os.path.join(analysis_dir, "crested", "enhancer_tables", "striatal_tools_supp_table1.xlsx"))
enhancer_table["enhancer ID"] = enhancer_table["enhancer ID"].str.strip()

## Coordinates
enhancer_table[["Chromosome", "Start", "End"]] = enhancer_table["Genome coordinates (cloned)"].apply(parse_coordinates)

## Sequence
enhancer_table["sequence"] = enhancer_table["Enhancer sequence"].str.lower()
enhancer_table["Sequence_CRESted"] = enhancer_table["sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
## Load in epi-conserved peaks
## -------------------------
epi_cons_markers = pd.read_csv(os.path.join(analysis_dir, "atac/conservation", "epi_conserved_marker_peaks.csv"))
epi_cons_markers[["Chromosome", "Start", "End"]] = epi_cons_markers["region"].apply(parse_coordinates)

## Run predictions
epi_cons_markers["sequence"] = [
   genome.fasta.fetch(chrom, start, end) 
   for chrom, start, end in tqdm(zip(epi_cons_markers["Chromosome"], epi_cons_markers["Start"], epi_cons_markers["End"]))
]
epi_cons_markers["Sequence_CRESted"] = epi_cons_markers["sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
## Load in the CRESted model
## -------------------------
model = tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False)
# model = tf.keras.models.load_model(os.path.join(model_dir, "human/99.keras"), compile=False)

## -------------------------
model_embedding = crested.tl.extract_layer_embeddings(
    input = epi_cons_markers["Sequence_CRESted"].tolist(),
    model = model,
    layer_name = "global_average_pooling1d",
)

## --------------------------
predictions = crested.tl.predict(
    input=epi_cons_markers["Sequence_CRESted"].tolist(),
    model=model,
    batch_size=64,
)

## Add to cellxgene as anndata
embedding_adata = ad.AnnData(
    X=predictions,
    var = pd.DataFrame(index=adata.obs_names),
    obsm={"X_embedding": model_embedding},
)
embedding_adata.obs["celltype_prediction_gini"] = calc_gini(embedding_adata.X).max(axis=1)
embedding_adata.obs["celltype_prediction"] = embedding_adata.var_names[embedding_adata.X.argmax(axis=1)]
embedding_adata.obs["neighborhood"] = embedding_adata.obs["celltype_prediction"].map(dict(zip(cluster_meta.Group, cluster_meta.Neighborhood)))

## -------------------------
## Plot UMAP
## -------------------------
plot_adata = embedding_adata[embedding_adata.obs["neighborhood"] != "Nonneuron"].copy()

sc.pp.neighbors(plot_adata, use_rep='X_embedding', n_neighbors=15)
sc.tl.umap(plot_adata, min_dist=0.3, random_state=42)

## Add to cellxgene as anndata
plot_adata.write_h5ad(os.path.join(cxg_dir, "CRESted", "basal-ganglia", "epi_conserved_peaks_crested_embedding.h5ad"))







## --- Build DataFrame ---
# embedding_df = pd.DataFrame(embedding, columns=[f"{method.upper()}1", f"{method.upper()}2"])
# embedding_df["celltype"] = embedding_adata["celltype_prediction"]
# embedding_df["gini_score"] = embedding_adata["celltype_prediction_gini"]

# ## --- Plotting ---
# plt.figure(figsize=(10, 10), dpi=900)
# for celltype in list(color_map.keys())[::-1]:
#     print(f"Plotting cell type: {celltype}")
#     if celltype in embedding_df["celltype"].values:
#         subset = embedding_df[embedding_df["celltype"] == celltype]
#         plt.scatter(
#             subset[f"{method.upper()}1"],
#             subset[f"{method.upper()}2"],
#             color=color_map[celltype],
#             alpha=0.6,
#             rasterized=True,
#             s=20,
#             label=celltype,
#         )
#         # Compute centroid for text label
#         x_center = subset[f"{method.upper()}1"].mean()
#         y_center = subset[f"{method.upper()}2"].mean()
#         # Add label at centroid
#         plt.text(
#             x_center,
#             y_center,
#             celltype,
#             fontsize=8,
#             ha='center',
#             va='center',
#             weight='bold',
#             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.3')
#         )

# ## Final plot
# plt.title(f"{method.upper()} Visualization of CRESTed Model Embeddings")
# plt.xlabel(f"{method.upper()}1")
# plt.ylabel(f"{method.upper()}2")
# plt.tight_layout()
# plt.savefig('/home/nelson.johansen/seq_umap_CERP_BG.png', bbox_inches='tight')

## Add to cellxgene as anndata
plot_adata.write_h5ad(os.path.join(cxg_dir, "CRESted", "basal-ganglia", "epi_conserved_peaks_crested_embedding.h5ad"))