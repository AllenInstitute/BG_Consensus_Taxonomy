import os
import sys
import crested
import keras
import glob
import re
import pandas as pd
import pyranges as pr
import anndata as ad
import tensorflow as tf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse
tf.config.list_physical_devices('GPU')

## -------------------------
## Parameters
species_name = "macaque" # sys.argv[1] ## species_name = "mouse"
run_name = f"{species_name}_basalganglia_hmba_pre-print" ## Specifies a unique directory for this model

## Paths
work_dir = "./aibs-octo-dnaseq-modeling/basal-ganglia/"
model_dir = os.path.join(work_dir, "models", "crested", species_name)

## -------------------------
if species_name in ["human", "mouse", "macaque", "marmoset"]:
    ## Zemke et al. 2023 data references
    species_to_reference = {
        "human": "hg38",
        "mouse": "mm10",
        "macaque": "rheMac10",
        "marmoset": "calJac4"
    }
    ## -------------------------
    reference = f"./aibs-octo-dnaseq-modeling/genomes/{species_name}"
    if species_name == "macaque":
        ## Macaque specific reference
        reference = os.path.join(reference, "ncbi")
        ##
        chrom_alias = pd.read_table(os.path.join(reference, "rheMac10.chromAlias.txt"), delimiter="\t")
        chrom_alias.columns = ["seqName", "aliasName", "UCSC"]
        ## Gather chr1-20, chrX,Y
        chromosomes = chrom_alias.loc[chrom_alias.UCSC.str.startswith("NC_")]
        chromosomes = chromosomes.loc[~chromosomes.seqName.str.contains("chrM")]
    ## Chromosome sizes
    chr_sizes = os.path.abspath(os.path.join(reference, f"{species_to_reference[species_name]}.chrom.sizes"))
    ## Fasta and GTF paths
    fasta_path = os.path.abspath(os.path.join(reference, f"fasta/genome.fa"))
    gtf_path = os.path.abspath(os.path.join(reference, f"genes/genes.gtf.gz"))
    ## Initalize genome object
    genome = crested.Genome(
        fasta=fasta_path,
        annotation=gtf_path,
    )

## Setup model directory
if ~os.path.exists(os.path.join(work_dir, "models", f"crested/{species_name}/")):
    os.makedirs(os.path.join(work_dir, "models", f"crested/{species_name}/"), exist_ok=True)

## Load data
adata = ad.read_h5ad(os.path.join(work_dir, "data", species_name, "crested_adata", f"{run_name}_finetune_crested.h5ad"))

##
datamodule = crested.tl.data.AnnDataModule(
    adata,
    genome=genome,
    batch_size=64,  ## A100 setting
    max_stochastic_shift=3, 
    always_reverse_complement=True,  # default True. Will double the effective size of the training dataset.
    chromsizes_file=chr_sizes,
)

##
model_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/macaque/99.keras"
model_architecture = keras.models.load_model(
    model_path, compile=False
)

## 
optimizer = keras.optimizers.Adam(learning_rate=1e-4)
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100, multiplier=1)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]
config = crested.tl.TaskConfig(optimizer, loss, metrics)

##
os.chdir(os.path.abspath(model_dir)) ## Kind of annoying that CREsted operates in a specific directory.
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=config,
    project_name="macaque_basal_ganglia_finetune",  ## Main directory for the project (species)
    run_name=run_name,  ## Use to distinguish different runs
    logger="wandb",  ## or 'wandb', 'tensorboard'
)

## train the model
trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=10,
    model_checkpointing_best_only=False
)

## ------------------------------
## Evaluate the model

## Run predictions
seqs = [
   genome.fasta.fetch(chrom, start, end) 
   for chrom, start, end in tqdm(zip(adata.var['chr'], adata.var['start'], adata.var['end']))
]

##
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/macaque/finetune"
model = keras.models.load_model(os.path.join(model_dir, "11.keras"), compile=False)

predictions = []
for batch in tqdm(range(0, len(seqs), 100000)):
    predictions.append(crested.tl.predict(seqs[batch:batch + 100000], model, batch_size=512, verbose=2))

predictions_ = np.concatenate(predictions, axis=0)
adata.layers["finetune"] = predictions_.T  # adata expects (C, N) instead of (N, C)

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