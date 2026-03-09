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
import scanpy as sc
from Bio import SeqIO

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
cactus_maf_dir = os.path.join(analysis_dir, "crested", "evo_analysis", "maf")

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
## Read in and process maf fasta
## -------------------------
maf_fasta_files = glob.glob(os.path.join(cactus_maf_dir, "*.fasta"))

## File path
fasta_file = maf_fasta_files[0]

## Dictionary to store species ID and sequence
sequences = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    species_id = record.id.split("/")[0]  ## get species id
    sequence = pad_sequence(str(record.seq).replace("-", "N"), target_width=2114) ## replace gaps with N's       
    sequences[species_id] = sequence

## -------------------------
## Load in the CRESted model
## -------------------------
# model = tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False)
model = tf.keras.models.load_model(os.path.join(model_dir, "human/99.keras"), compile=False)

## -------------------------
# model_embedding = crested.tl.extract_layer_embeddings(
#     input =adata_top.var["sequence"].tolist(),
#     model = model,
#     layer_name = "global_average_pooling1d",
# )

## --------------------------
predictions = crested.tl.predict(
    input=list(sequences.values()),
    model=model,
    batch_size=64,
)

## Add to cellxgene as anndata
embedding_adata = ad.AnnData(
    X=predictions.T,
    obs = pd.DataFrame(index=adata.obs_names),
    var = pd.DataFrame(index=sequences.keys()),
    # obsm={"X_embedding": model_embedding},
    # layers={"accessibility": adata.X.T},
)

## Group order
group_order = cluster_meta["Group"][cluster_meta["Group"].isin(embedding_adata.obs_names)]

## Species order
species_order = cactus_adata.var.sort_values('primates_shorthand').index
species_order = species_order[species_order.isin(embedding_adata.var_names)]

## Plot
plt.figure(figsize=(10,8))
sns.heatmap(embedding_adata[group_order, species_order].X, 
            annot=False, 
            fmt=None,
            cmap="Reds",
            cbar=True,
            xticklabels=False, 
            yticklabels=group_order)
plt.xlabel("Species (Alignment)")
plt.ylabel("Group (Cell Type)")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "crested", "evo_analysis", "figures", f"human_maf_alignment_predictions.png"), dpi=300)

##
ax = sc.pl.heatmap(embedding_adata, groupby="Group", cmap="viridis", dendrogram=True)