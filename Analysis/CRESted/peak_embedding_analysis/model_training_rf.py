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

from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import precision_recall_fscore_support as scores

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
species = "macaque"
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
top_k = 100  # Here we take the top 500 most specific regions
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
results = {}
for model_type in ["pretrained", "finetune"]:  #"pretrained",
    if model_type == "finetune":
        model_files = glob.glob(os.path.join(model_dir, f"{species}/{model_type}/*.keras"))
    else:
        model_files = glob.glob(os.path.join(model_dir, f"{species}/*.keras"))
    ## Sort model files from 0->max epoch, model files are named #.keras
    model_files = sorted(model_files, key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    ## Run through models and do random forest classification on embeddings for top 100 peaks
    for model_path in model_files:
        print(f"Processing model: {model_path}")
        model = tf.keras.models.load_model(model_path, compile=False)
        ## -------------------------
        model_embedding = crested.tl.extract_layer_embeddings(
            input=adata_top.var["sequence"].tolist(),
            model=model,
            layer_name="global_average_pooling1d",
        )
        ## Add to cellxgene as anndata
        embedding_adata = ad.AnnData(
            X=model_embedding,
            obs=adata_top.var.copy(),
            obsm={"X_embedding": model_embedding},
        )
        embedding_adata.obs_names_make_unique()
        embedding_adata.obs["celltype_prediction_gini"] = calc_gini(embedding_adata.X).max(axis=1)
        embedding_adata.obs["celltype_prediction"] = embedding_adata.var_names[embedding_adata.X.argmax(axis=1)]
        ##
        for tax_level in ["Neighborhood", "Class", "Subclass"]:
            embedding_adata.obs[tax_level.lower()] = embedding_adata.obs["celltype"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))
            embedding_adata.obs[f"{tax_level.lower()}_pred"] = embedding_adata.obs["celltype_prediction"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))
        ## -------------------------
        ## Accuracy of cell type prediction
        ## -------------------------
        ## Initialize Random Forest
        clf = RandomForestClassifier(
            n_estimators=200,    # number of trees (can tune)
            max_depth=5,     # let trees grow fully (or set a limit)
            random_state=42,
            n_jobs=-1           # use all CPU cores for speed
        )
        ## Train on all data (just for reporting)
        clf.fit(embedding_adata.X, embedding_adata.obs["celltype"])
        embedding_adata.obs["celltype_pred_rf"] = clf.predict(embedding_adata.X)
        ## Return the per cell type F1 score
        p, r, f1, s = scores(embedding_adata.obs["celltype"], embedding_adata.obs["celltype_pred_rf"],
                        labels=embedding_adata.obs["celltype"].astype("category").cat.categories.to_list(),
                        average=None,
                        sample_weight=None)
        result = pd.DataFrame({'label': embedding_adata.obs["celltype"].astype("category").cat.categories,
                            'precision': p,
                            'recall': r,
                            'f1': f1,
                            'support': s,
                            'epoch': int(re.search(r'(\d+)', os.path.basename(model_path)).group(1)),
                            'model_type': model_type})
        print(result.epoch[0].astype(str) + f"_{model_type}:" + f1.mean().astype(str))
        results[result.epoch[0].astype(str) + f"_{model_type}"] = result

## -------------------------
## Store results and plot
results_df = pd.concat(results.values(), axis=0)
results_df["step"] = results_df["epoch"].astype(str) + "_" + results_df["model_type"]

## Compute mean F1 per model
mean_f1 = results_df.groupby("step")["f1"].mean().reset_index()

## Add back into the dataframe
results_df = results_df.merge(mean_f1, on="step", suffixes=("", "_mean"))
results_df_per_step = results_df.drop_duplicates(subset=["step"])[["step", "f1_mean"]].reset_index(drop=True)

## -------------------------
## Plot as a lineplot with dots, highlight the pretrain and finetune steps
plt.figure(figsize=(4,3))
sns.lineplot(data=results_df_per_step, x="step", y="f1_mean", marker="o", markersize=2.5, markeredgewidth=0.25, markeredgecolor="black", linewidth=0.5)
plt.xticks(rotation=45, ha='right', fontsize=4)
plt.xlabel("Model Checkpoint (Epoch_ModelType)")
plt.ylabel("Mean F1 Score Across Cell Types")
plt.title(f"Random Forest Classification Performance on {top_k} Most Specific Regions - {species.capitalize()}")
plt.grid(
    linestyle='-',  
    linewidth=0.25, 
    alpha=0.4
)
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, f"crested/figures_qc/{species}_rf_classification_performance_top{top_k}_regions.pdf"), dpi=900)
plt.close()

## -------------------------
## Plot confusion matrix for best model
cm = confusion_matrix(embedding_adata.obs["celltype"], embedding_adata.obs["celltype_pred_rf"], labels=clf.classes_)

## Reorder cm according to cluster_meta
ordered_classes = cluster_meta["Group"][cluster_meta["Group"].isin(clf.classes_)].values
ordered_indices = [list(clf.classes_).index(cls) for cls in ordered_classes]
cm = cm[np.ix_(ordered_indices, ordered_indices)]

plt.figure(figsize=(8, 8))
sns.heatmap(
    cm,
    cmap='Blues',
    xticklabels=ordered_classes,
    yticklabels=ordered_classes,
    cbar=False
)
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Confusion Matrix — Random Forest")
plt.xticks(rotation=45, ha='right', fontsize=4)
plt.yticks(rotation=0, fontsize=4)
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, f"crested/figures_qc/{species}_rf_classification_confusion_matrix_top{top_k}_regions.pdf"), dpi=900)