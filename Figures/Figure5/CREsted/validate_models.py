import os, re, sys, glob
import anndata as ad
import pandas as pd
import keras
import crested
import pysam
import numpy as np
from tqdm import tqdm

import tensorflow as tf

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print("Memory growth enabled for GPUs")
    except RuntimeError as e:
        print("Failed to set memory growth:", e)

## --------------------------
## Load model
species_name = "human"
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling"

## -------------------------
## Setup genome
reference = os.path.join(work_dir, f"genomes/{species_name}")
if species_name == "macaque":
    reference = os.path.join(reference, "ncbi")

## Chromosome sizes
if os.path.exists(os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))):
    chr_sizes = os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))
else:
    chr_sizes = None

## Fasta and GTF paths
fasta_path = os.path.abspath(os.path.join(reference, f"fasta/genome.fa"))
fasta = pysam.FastaFile(fasta_path)

## ---------------------------
## Load the pretrained CREsted model
model_name = "CrestedModelAPI--dilated-cnn-human-basalganglia--25-06-05-13-57-57"

## Check if model is finetuned or not.
model_dir = os.path.join(work_dir, "basal-ganglia", "models", f"crested/evogen-dnaseq-modeling/basal_ganglia/{model_name}")
if os.path.exists(model_dir + "_finetune"):
    model_dir = model_dir + "_finetune/checkpoints/"
else:
    model_dir = model_dir + "/checkpoints/"

## Find the model with highest training number this isn't IDEAL as the best model may not be the last one
best_model = max(glob.glob(os.path.join(model_dir, "*.keras")), key=lambda x: int(re.search(r'(\d+)\.keras$', x).group(1)))
print("Loading best model: ", best_model)

## Initalize model
model = keras.models.load_model(best_model, compile=False)

## ---------------------------
## Gather data model was trained on, really only needed to determine class (celltype) names
adata = ad.read_h5ad(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata", f"{species_name}_basalganglia_hmba_pre-print_crested.h5ad"))

## Run predictions
seqs = [
   genome.fasta.fetch(chrom, start, end) 
   for chrom, start, end in tqdm(zip(adata.var['chr'], adata.var['start'], adata.var['end']))
]

predictions = []
for batch in tqdm(range(0, len(seqs), 100000)):
    predictions.append(crested.tl.predict(seqs[batch:batch + 100000], model, batch_size=512, verbose=2))

predictions_ = np.concatenate(predictions, axis=0)
adata.layers["base"] = predictions_.T  # adata expects (C, N) instead of (N, C)

## ---------------------------
## Setup output directories
figure_dir = os.path.join("/home/nelson.johansen/", "striatal_tools", "figures", f"{species_name}")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

##
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables"
consensus_anno = pd.read_csv(os.path.join(anno_dir, "consensus_anno_pre-print.csv"), index_col=0)

##
consensus_anno["Group"] = consensus_anno["Group"].str.replace(" ", "_").unique()
consensus_anno = consensus_anno.loc[consensus_anno["Group"].isin(adata.obs_names.unique()), :]

## Re order data object based on consensus annotation
adata = adata[consensus_anno["Group"].values, :]

## ---------------------------
## Run enhancer sequences through model
crested.pl.heatmap.correlations_predictions(
    adata,
    split="test",
    title="Correlations between Groundtruths and Predictions",
    x_label_rotation=90,
    width=15,
    height=15,
    log_transform=True,
    reorder=False,
    vmax=1,
    vmin=-0.15,
    save_path=os.path.join(figure_dir, f"{species_name}_basalganglia_crested_correlations_predictions.png"),
)

crested.pl.heatmap.correlations_self(
    adata,
    title="Self Correlation Heatmap",
    x_label_rotation=90,
    width=15,
    height=15,
    log_transform=True,
    reorder=False,
    vmax=1,
    vmin=-0.15,
    save_path=os.path.join(figure_dir, f"{species_name}_basalganglia_crested_self_correlations.png"),
)