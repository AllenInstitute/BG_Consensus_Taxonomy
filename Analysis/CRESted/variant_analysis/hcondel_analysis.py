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
variant_lib = pd.read_excel(os.path.join(variant_dir, "hCONDELs_10032_2023Science.xlsx"), sheet_name="hCONDEL metadata")

## Genome 
fasta_ = pysam.FastaFile(os.path.join(work_dir, "genomes", "human", "fasta", "hg38.fa"))

## Filter for enhancers
variant_lib["Sequence"] = None
variant_lib["Sequence"] = variant_lib.apply(
    lambda row: fasta_.fetch(row["hg38_cons_chr"], row["hg38_cons_start_pos"], row["hg38_cons_end_pos"]),
    axis=1
)

variant_lib["del_Sequence"] = None

## Create deletion sequences by masking out deleted region with N's
for row in range(len(variant_lib)):
    ref_seq = list(str(variant_lib["Sequence"].values[row]))
    ## Find overlap between subsequence and deletion
    overlap_start = max(variant_lib.loc[row, "hg38_cons_start_pos"], variant_lib.loc[row, "hg38_del_start_pos"]) - 1 
    overlap_end   = min(variant_lib.loc[row, "hg38_cons_end_pos"], variant_lib.loc[row, "hg38_del_end_pos"]) - 1
    ## Convert to local coordinates
    local_start = max(0, overlap_start - variant_lib.loc[row, "hg38_cons_start_pos"])
    local_end   = min(len(ref_seq) - 1, overlap_end - variant_lib.loc[row, "hg38_cons_start_pos"])
    for i in range(local_start, local_end + 1):
        ref_seq[i] = "N"
    ##
    variant_lib.loc[row, "del_Sequence"] = "".join(ref_seq)

## --------------------------
## Load model
## --------------------------

## ---------------------------
## Load the pretrained CREsted model

variant_lib["Sequence_CRESted"] = variant_lib["Sequence"].apply(lambda x: pad_sequence(x, 2114))
variant_lib["del_Sequence_CRESted"] = variant_lib["del_Sequence"].apply(lambda x: pad_sequence(x, 2114))

## -------------------------
model = tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/11.keras"), compile=False)
adata = ad.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/human/crested_adata/human_basalganglia_hmba_pre-print_finetune_crested.h5ad")

## Predict for eqch consensus sequence
predictions = crested.tl.predict(
    input=variant_lib["Sequence_CRESted"].tolist(),
    model=model,
    batch_size=16,
)
predictions_df = pd.DataFrame(predictions, 
                                index=variant_lib["hCONDEL_ID"], 
                                columns=adata.obs_names)


## Predict for each deletion sequence
predictions_del = crested.tl.predict(
    input=variant_lib["del_Sequence_CRESted"].tolist(),
    model=model,
    batch_size=16,
)
predictions_del_df = pd.DataFrame(predictions_del, 
                                index=variant_lib["hCONDEL_ID"], 
                                columns=adata.obs_names)

##
difference_df = predictions_df - predictions_del_df
difference_df.to_csv(os.path.join(figure_dir, "hCONDEL_CRESted_predictions_difference.csv"))

##
fold_change_df = np.log2((predictions_df + 1e-8) / (predictions_del_df + 1e-8))
fold_change_df.to_csv(os.path.join(figure_dir, "hCONDEL_CRESted_predictions_log2FC.csv"))



# Plot
plt.figure(figsize=(10,8))
predictions_df.T[0].plot(kind="bar")

plt.xlabel("Category")
plt.ylabel("Value")
plt.title("Barplot from 1-column DataFrame")
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, "CRESted_predictions_example.png"), dpi=300)
plt.close()