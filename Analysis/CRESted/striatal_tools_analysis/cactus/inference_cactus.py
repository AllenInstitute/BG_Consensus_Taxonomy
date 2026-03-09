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

import pickle
from Bio import SeqIO

# === TensorFlow Setup ===
tf.config.list_physical_devices('GPU')

## Function to pad sequences equally on both sides
def pad_sequence(sequence, target_width):
    ## If sequence is to long cut out the middle
    if len(sequence) >= target_width:
        centers = len(sequence) // 2
        start = int(centers - (target_width / 2))
        end = int(centers + (target_width / 2))
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

## Location of maf files
cactus_maf_dir = os.path.join(analysis_dir, "crested", "maf")

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
cactus_adata = ad.read_h5ad(os.path.join(analysis_dir, "cactus", "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"), backed="r")

## -------------------------
## Load in striatal tools table
## -------------------------
striatal_table = pd.read_excel(os.path.join(analysis_dir, "crested", "enhancer_tables", "striatal_tools_supp_table1.xlsx"))
striatal_table["enhancer ID"] = striatal_table["enhancer ID"].str.strip()

## Coordinates
striatal_table[["Chromosome", "Start", "End"]] = striatal_table["Genome coordinates (cloned)"].apply(parse_coordinates)
## Sequence
striatal_table["sequence"] = striatal_table["Enhancer sequence"].str.lower()
striatal_table["Sequence_CRESted"] = striatal_table["sequence"].apply(lambda x: pad_sequence(x, 2114))

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
## Read in and process maf fasta
## -------------------------
maf_fasta_files = glob.glob(os.path.join(cactus_maf_dir, "*_stitched_and_aligned.fasta"))

maf_enhancer_sequences = {}

for fasta_file in maf_fasta_files:
    enhancer_id = re.search(r"([^/]+)_stitched_and_aligned\.fasta$", fasta_file).group(1)
    ## Dictionary to store species ID and sequence
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        species_id = record.id.split("/")[0]  ## get species id
        sequence = pad_sequence(str(record.seq).replace("-", "N"), target_width=2114) ## replace gaps with N's       
        sequences[species_id] = sequence
    maf_enhancer_sequences[enhancer_id] = sequences

## Convert to a 3 column dataframe with enhancer_id, species, sequence
enhancer_table = []
for enhancer_id, seq_dict in maf_enhancer_sequences.items():
    for species_id, sequence in seq_dict.items():
        enhancer_table.append({
            "enhancer_id": enhancer_id,
            "species": species_id,
            "sequence": sequence
        })

enhancer_table = pd.DataFrame(enhancer_table)

enhancer_table["Sequence_CRESted"] = enhancer_table["sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
## Load in the CRESted model
## -------------------------
model = tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False)

## -------------------------
model_embedding = crested.tl.extract_layer_embeddings(
    input = enhancer_table["Sequence_CRESted"].tolist(),
    model = model,
    layer_name = "global_average_pooling1d",
)

## --------------------------
predictions = crested.tl.predict(
    input=enhancer_table["Sequence_CRESted"].tolist(),
    model=model,
    batch_size=64,
)

## Add to cellxgene as anndata with empty X matrix filled with 0
pred_maf_adata = ad.AnnData(
    X=np.zeros((len(enhancer_table["enhancer_id"].unique()), len(cactus_adata.var_names))),
    obs = pd.DataFrame(index=enhancer_table["enhancer_id"].unique()), #enhancer_table,
    var = pd.DataFrame(index=cactus_adata.var_names)
)

## Iterate through each enhancer and species to fill in the X matrix 
for itr, row in tqdm(enhancer_table.iterrows()):
    enhancer_id = row["enhancer_id"]
    species_id = row["species"]
    prediction = predictions[itr]
    ## 
    celltype = striatal_table.loc[striatal_table["enhancer ID"] == enhancer_id, "Cell type"].values[0]
    if celltype == "Pan-MSN":
        celltype_idx = adata.obs_names.str.contains("MSN")
    elif celltype == "D1 MSN":
        celltype_idx = adata.obs_names.str.contains("D1") & adata.obs_names.str.contains("MSN")
    elif celltype == "D2 MSN":
        celltype_idx = adata.obs_names.str.contains("D2") & adata.obs_names.str.contains("MSN")
    elif celltype == "Cholinergic":
        celltype_idx = adata.obs_names.str.contains("Cholinergic")
    elif celltype == "Sst-Chodl":
        celltype_idx = adata.obs_names.str.contains("SST-CHODL")
    elif celltype == "Pvalb-Pthlh":
        celltype_idx = adata.obs_names.str.contains("PTHLH-PVALB")
    ## Store embedding and prediction in the corresponding obs
    pred_maf_adata[enhancer_id, species_id].X = prediction[celltype_idx].mean()

## Write out
with open(os.path.join(analysis_dir, "crested", "data", "striatal_tools_enhancer_maf_predictions.pkl"), "wb") as f:
    pickle.dump(predictions, f)

##
pred_maf_adata.write_h5ad(os.path.join(analysis_dir, "crested", "data", "striatal_tools_enhancer_maf_predictions_v2.h5ad"))

