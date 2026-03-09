import numpy as np 
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import keras
import crested
import anndata as ad
import os
import tensorflow as tf
import pysam

## -------------------------
## Helpful paths
## -------------------------
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
## Model related
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
## Load models and associated data
## -----------------------
models = {
    "human": tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False),
    "macaque": tf.keras.models.load_model(os.path.join(model_dir, "macaque/finetune/11.keras"), compile=False)
}

adatas = {
    "human": ad.read_h5ad(os.path.join(data_dir, "human/crested_adata/human_basalganglia_hmba_pre-print_finetune_crested.h5ad")),
    "macaque": ad.read_h5ad(os.path.join(data_dir, "macaque/crested_adata/macaque_basalganglia_hmba_pre-print_finetune_crested.h5ad"))
}

## Genomes
genome_human = "/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/fasta/genome.fa"
genome_macaque  = "/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/fasta/genome.fa"

## -------------------------
## Run predictions
## -------------------------
predictions = crested.tl.predict(adatas["human"], models["human"], batch_size=128, genome=genome_human)
adatas["human"].layers["finetuned"] = predictions.T

predictions = crested.tl.predict(adatas["macaque"], models["macaque"], batch_size=128, genome=genome_macaque)
adatas["macaque"].layers["finetuned"] = predictions.T

## -------------------------
## Save adatas with predictions
## -------------------------
adatas["human"].write_h5ad(os.path.join(analysis_dir, f"crested/figures_qc/data/human_basalganglia_hmba_pre-print_finetune_crested_model_qc.h5ad"))
adatas["macaque"].write_h5ad(os.path.join(analysis_dir, f"crested/figures_qc/data/macaque_basalganglia_hmba_pre-print_finetune_crested_model_qc.h5ad"))

## Reorder adata to match taxonomy
for species in adatas.keys():
    adatas[species] = adatas[species][cluster_meta["Group"][cluster_meta["Group"].isin(adatas[species].obs_names)].values]

## -------------------------
## Plot model predictions vs targets
## -------------------------
for species, adata in adatas.items():
    with plt.rc_context({'figure.figsize': (12,12), 'figure.dpi': 900}):
        crested.pl.scatter.class_density(
            adata,
            class_name=None,
            model_names=["finetuned"],
            split="test",
            log_transform=True,
            density_indication=True,
            width=5,
            height=5,
            title=species + ' model predictions vs targets',
            save_path=os.path.join(analysis_dir, f"crested/figures_qc/{species}_finetune_model_predictions_vs_targets.png")
        )

## -------------------------
for species, adata in adatas.items():
    with plt.rc_context({'figure.figsize': (16,18), 'figure.dpi': 900}):
        crested.pl.heatmap.correlations_self(
            adata, 
            log_transform=True, 
            vmin=0.4,
            vmax=0.85,
            x_label_fontsize=6,
            y_label_fontsize=6,
            x_label_rotation=90,
            title="Self correlations heatmap",
            save_path=os.path.join(analysis_dir, f"crested/figures_qc/{species}_finetune_model_self_correlation_heatmap.png")
        )
        crested.pl.heatmap.correlations_predictions(
            adata,
            model_names=["finetuned"],
            split="test",
            log_transform=True,
            vmin=0.4,
            vmax=0.85,
            x_label_fontsize=6,
            y_label_fontsize=6,
            x_label_rotation=90,
            title="Correlations: Predictions vs Ground Truth",
            save_path=os.path.join(analysis_dir, f"crested/figures_qc/{species}_finetune_model_correlations_predictions.png")
        )