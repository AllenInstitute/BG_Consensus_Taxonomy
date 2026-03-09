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
species = "human"
adata = ad.read_h5ad(os.path.join(data_dir, f"{species}/crested_adata/{species}_basalganglia_hmba_pre-print_crested.h5ad"))

## -------------------------
## Genome
## -------------------------
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
## Gather top 500
## -------------------------

## Save the original data
adata_top = adata.copy()

## Sort and filter to top 500 most specific regions
top_k = 500  # Here we take the top 500 most specific regions
crested.pp.sort_and_filter_regions_on_specificity(
    adata_top, top_k=top_k, method="proportion"
)
adata_top.var.rename(columns={"Class name": "celltype"}, inplace=True)

## Gather sequences for top 500 regions which are the right length for CRESted already!
adata_top.var["sequence"] = [
   genome.fasta.fetch(chrom, start, end) 
   for chrom, start, end in tqdm(zip(adata_top.var["chr"], adata_top.var["start"], adata_top.var["end"]))
]

## -------------------------
## Load in the CRESted model
## -------------------------
model_type = "finetune"
if model_type == "finetune":
    model = tf.keras.models.load_model(os.path.join(model_dir, f"{species}/finetune/12.keras"), compile=False)
else:
    model = tf.keras.models.load_model(os.path.join(model_dir, f"{species}/99.keras"), compile=False)

## -------------------------
model_embedding = crested.tl.extract_layer_embeddings(
    input =adata_top.var["sequence"].tolist(),
    model = model,
    layer_name = "global_average_pooling1d",
)

## --------------------------
predictions = crested.tl.predict(
    input=adata_top.var["sequence"].tolist(),
    model=model,
    batch_size=64,
)

## Add to cellxgene as anndata
embedding_adata = ad.AnnData(
    X=predictions,
    obs = adata_top.var.copy(),
    var = pd.DataFrame(index=adata_top.obs_names),
    obsm={"X_embedding": model_embedding},
    layers={"accessibility": adata_top.X.T},
)
embedding_adata.obs_names_make_unique()
embedding_adata.obs["celltype_prediction_gini"] = calc_gini(embedding_adata.X).max(axis=1)
embedding_adata.obs["celltype_prediction"] = embedding_adata.var_names[embedding_adata.X.argmax(axis=1)]

for tax_level in ["Neighborhood", "Class", "Subclass"]:
    embedding_adata.obs[tax_level.lower()] = embedding_adata.obs["celltype"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))
    embedding_adata.obs[f"{tax_level.lower()}_pred"] = embedding_adata.obs["celltype_prediction"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))

## Compute correlation between gini_score and celltype_prediction_gini
# from scipy.stats import pearsonr, spearmanr
# corr_pearson, pval_pearson = pearsonr(embedding_adata.obs["gini_score"], embedding_adata.obs["celltype_prediction_gini"])
# corr_spearman, pval_spearman = spearmanr(embedding_adata.obs["gini_score"], embedding_adata.obs["celltype_prediction_gini"])
# print(f"Pearson correlation: {corr_pearson:.4f} (p={pval_pearson:.4e})")
# print(f"Spearman correlation: {corr_spearman:.4f} (p={pval_spearman:.4e})")

## -------------------------
## Plot UMAP
## -------------------------
plot_adata = embedding_adata[~((embedding_adata.obs["neighborhood"] == "Nonneuron") |
                             (embedding_adata.obs["neighborhood_pred"] == "Nonneuron"))].copy()
# plot_adata = plot_adata[(plot_adata.obs["proportion_score"] > 0.12) & (plot_adata.obs["celltype_prediction_gini"] > 0.4)].copy()

sc.pp.neighbors(plot_adata, use_rep='X_embedding', n_neighbors=15)
sc.tl.umap(plot_adata, min_dist=0.3, random_state=42)

## Add to cellxgene as anndata
# plot_adata.write_h5ad(os.path.join(cxg_dir, "CRESted", "basal-ganglia", f"top_peaks_100_{species}_crested_{model_type}_embedding.h5ad"))

## -------------------------
## UMAP figure
## -------------------------
## Setup colors for plotting
colors = plot_adata.obs["celltype"].map(color_map).values.astype(str)
colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values

## Plot UMAP colored by cell type
plt.figure(figsize=(12,12))
plt.scatter(plot_adata.obsm["X_umap"][:,0], 
            plot_adata.obsm["X_umap"][:,1], 
            c=colors,
            s=15, 
            alpha=0.8)
## Add color legend
plt.title("UMAP of CRESted embedding of top 100 Human peaks")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "crested", "figures", f"{species}_{model_type}_crested_umap_top100.png"), dpi=900)  

## -------------------------
## Accuracy of cell type prediction
## -------------------------

from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

## Initialize Random Forest
clf = RandomForestClassifier(
    n_estimators=200,    # number of trees (can tune)
    max_depth=None,     # let trees grow fully (or set a limit)
    random_state=42,
    n_jobs=-1           # use all CPU cores for speed
)

## Train on all data (just for reporting)
clf.fit(embedding_adata.X, embedding_adata.obs["celltype"])
embedding_adata.obs["celltype_pred_rf"] = clf.predict(embedding_adata.X)

## Confusion matrix
cm = confusion_matrix(embedding_adata.obs["celltype"], 
                      embedding_adata.obs["celltype_pred_rf"], 
                      labels=cluster_meta["Group"].tolist())

## -- Plot
plt.figure(figsize=(10,8))
sns.heatmap(cm/100, 
            annot=False, 
            fmt=None,
            cmap="Blues",
            cbar=True,
            xticklabels=cluster_meta["Group"].tolist(), 
            yticklabels=cluster_meta["Group"].tolist())
plt.xlabel("Predicted")
plt.ylabel("True")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "crested", "figures", f"{species}_{model_type}_top500_celltype_confusion_matrix.png"), dpi=300)

## Save the embedding adata
embedding_adata.write_h5ad(os.path.join(analysis_dir, "crested", "data", f"{species}_top500_crested_embedding.h5ad"))

## -------------------------
## Plot contribution scores for the top 5 peaks per cell type
## -------------------------


for celltype in plot_adata.obs.celltype.unique():
    print(f"Processing cell type: {celltype}")
    ## Get top 5 peaks for this cell type
    sub_adata = plot_adata[plot_adata.obs.celltype == celltype]
    peak = sub_adata.obs.sort_values(by="celltype_prediction_gini", ascending=False).index[0]
    ##
    seq = plot_adata.obs.loc[peak, "sequence"]
    cl_idx = list(plot_adata.var_names).index(celltype)
    ##
    regions_of_interest=[seq]
    scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
        regions_of_interest,
        target_idx=[cl_idx],
        model=model,
        method='integrated_grad'
    )
    ##
    crested.pl.patterns.contribution_scores(
        scores,
        one_hot_encoded_sequences,
        sequence_labels=[""],
        class_labels=[np.array(plot_adata.obs_names)[cl_idx]],
        zoom_n_bases=300,
        title=f"{celltype} - {peak}",
        height=5,
        save_path=os.path.join(analysis_dir, "crested", "figures", f"{species}_{model_type}_contribution_{celltype}_{peak}.pdf")
    )  # zoom in on the center 500bp
