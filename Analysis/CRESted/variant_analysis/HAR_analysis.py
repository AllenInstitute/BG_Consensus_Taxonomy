import os
import sys
import glob
import re
import anndata as ad
import pandas as pd
from pysam import FastaFile

import keras
import crested
import tensorflow as tf

## Function to pad sequences equally on both sides
def pad_sequence(sequence, target_width):
    ## If sequence is to long cut out the middle
    if len(sequence) >= target_width:
        center = len(sequence) // 2
        start = int(center - (target_width / 2))
        end = int(center + (target_width / 2))
        return sequence[start:end]
    ## Otherwise pad the sequence with N's
    total_padding = target_width - len(sequence)
    left_padding = total_padding // 2
    right_padding = total_padding - left_padding
    return ('N' * left_padding) + sequence + ('N' * right_padding)

##
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/figures/CRESted"
species_name = "human"

## Load enhancer library
variant_lib = pd.read_excel(os.path.join(variant_dir, "HARs_312_hg38_Science2023.xlsx"), sheet_name="zooHARs")

## Genome 
fasta_ = pysam.FastaFile(os.path.join(work_dir, "genomes", "human", "fasta", "hg38.fa"))

## Filter for enhancers
variant_lib["Sequence"] = None
variant_lib["Sequence"] = variant_lib.apply(
    lambda row: fasta_.fetch(row["chrom"], row["start"], row["end"]),
    axis=1
)

variant_lib["length"] = variant_lib["end"] - variant_lib["start"]

## --------------------------
## Load model
## --------------------------

## ---------------------------
## Load the pretrained CREsted model

variant_lib["Sequence_CRESted"] = variant_lib["Sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
species_name = "macaque"
model = tf.keras.models.load_model(os.path.join(model_dir, f"{species_name}/finetune/11.keras"), compile=False)
adata = ad.read_h5ad(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/{species_name}/crested_adata/{species_name}_basalganglia_hmba_pre-print_finetune_crested.h5ad")

## Predict for eqch consensus sequence
predictions = crested.tl.predict(
    input=["".join(random.choices(["A", "T", "C", "G"], k=2114))], #variant_lib["Sequence_CRESted"].tolist(),
    model=model,
    batch_size=16,
)
predictions_df = pd.DataFrame(predictions, 
                                #index=variant_lib["simple_name"], 
                                columns=adata.obs_names)
predictions_df.to_csv(os.path.join(figure_dir, "HAR_human_CRESted_predictions.csv"))

# Plot
plt.figure(figsize=(10,8))
predictions_df.T[0].plot(kind="bar")

plt.xlabel("Category")
plt.ylabel("Value")
plt.title("Barplot from 1-column DataFrame")
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, "CRESted_predictions_example_ATCG.png"), dpi=300)
plt.close()

## --------------------------
## Contribtion scores for HARs

scores, seqs_one_hot = crested.tl.contribution_scores(
    input=variant_lib["Sequence_CRESted"].tolist(),
    target_idx=None,
    model = model,
    method = "expected_integrated_grad"
)


def zoom_and_mean(scores, seqs_one_hot, zoom_levels, celltypes, ids):
    """
    scores: 4D array (batch, channels, length, features)
    zoom_levels: list of window sizes (int)
    
    returns: mean over zooms (same shape as scores, except length replaced by len(zoom_levels))
    """
    center = scores.shape[2] // 2
    zoomed_slices = []
    for idx, zoom in enumerate(zoom_levels):
        start_idx = center - zoom // 2
        end_idx = start_idx + zoom
        ## slice but keep 4D (length dimension shrinks)
        zoomed_scores = scores[idx, :, start_idx:end_idx, :]
        zoomed_mask = seqs_one_hot[idx, start_idx:end_idx, :]
        reference_scores = zoomed_scores * zoomed_mask[None,:]
        reference_scores = reference_scores.sum(axis=-1) ## Collapse last dimension which only has 1 value in the reference nucleotide, output: seq x cell_types x positions
        reference_scores = reference_scores.mean(axis=-1) ## Mean over the nucleotide positions and collapse, output: seq x cell_types
        ##
        zoomed_slices.append(reference_scores)  ## Add back a length dimension
    ## Turn zoomed_slices into a pandas dataframe
    zoomed_slices = pd.DataFrame(zoomed_slices, columns=celltypes, index=ids)
    return zoomed_slices

##
summary_scores = zoom_and_mean(scores, 
                                seqs_one_hot, 
                                variant_lib["length"],
                                adata.obs_names.tolist(), 
                                variant_lib["simple_name"].tolist())
summary_scores.to_csv(os.path.join(figure_dir, "HAR_human_CRESted_contribution_scores_summary.csv"))
# ##
# crested.pl.patterns.contribution_scores(
#     scores[0:1,:,:], 
#     seqs_one_hot[0:1,:,:], 
#     sequence_labels = variant_lib["simple_name"].tolist()[0:1], 
#     class_labels = adata.obs_names.tolist(),
#     zoom_n_bases=variant_lib["length"][0],
# )
# plt.savefig(os.path.join(figure_dir, "CRESted_contribution_scores_example.png"), dpi=300)