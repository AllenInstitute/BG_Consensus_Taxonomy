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
from tqdm import tqdm
tf.config.list_physical_devices('GPU')

import pyranges as pr
from pyfaidx import Fasta


## -------------------------
## Parameters
species_name = "marmoset" # sys.argv[1] ## species_name = "mouse"
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
    elif species_name == "marmoset":
        reference = os.path.join(reference, "hmba")
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
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia"
adata = ad.read_h5ad(os.path.join(work_dir, "data", species_name, "crested_adata", f"{run_name}_crested_10925.h5ad"))

##
os.chdir("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/marmoset")

##
datamodule = crested.tl.data.AnnDataModule(
    adata,
    genome=genome,
    batch_size=128,  ## A100 setting
    max_stochastic_shift=3, 
    always_reverse_complement=True,  # default True. Will double the effective size of the training dataset.
    chromsizes_file=chr_sizes,
)

## Clear the backend just encase a graph is still active
tf.keras.backend.clear_session()

## Model architecture
model_architecture = crested.tl.zoo.chrombpnet(
    seq_len=2114, 
    num_classes=len(list(adata.obs_names))
)

## Load the default configuration for training a peak regression model
config = crested.tl.default_configs(
    "peak_regression"
) 

##
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=config,
    project_name=species_name,  ## Main directory for the project (species)
    run_name=run_name,  ## Use to distinguish different runs
    logger="tensorboard",  ## or 'wandb', 'tensorboard'
)

## train the model
trainer.fit(epochs=100,
            model_checkpointing_best_only=False,
            early_stopping=False,
            learning_rate_reduce_patience=10)



## ------------------------------
## Evaluate the model

##
model_path = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/{species_name}/marmoset/10.keras"
# model_architecture = keras.models.load_model(
#     model_path, compile=False
# )

## Run predictions
# seqs = [
#    genome.fasta.fetch(chrom, start, end) 
#    for chrom, start, end in tqdm(zip(adata.var['chr'], adata.var['start'], adata.var['end']))
# ]


## Function to fetch sequences in batch
def fetch_sequences(chrom, starts, ends):
    ref = fasta[chrom]
    return [ref[start:end].seq for start, end in zip(starts, ends)]

## Load genome FASTA
fasta = Fasta(fasta_path)

## Convert BED to PyRanges
gr = pr.PyRanges(adata.var.rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'}))

# Add index column
gr = gr.df
gr["orig_idx"] = gr.index
gr = pr.PyRanges(gr)

seqs_per_idx = {}
for chrom, subdf in tqdm(gr.df.groupby("Chromosome")):
    if chrom not in fasta.keys():
        seqs_chr = ["N" * 500] * len(subdf)
    else:
        seqs_chr = fetch_sequences(chrom, subdf["Start"], subdf["End"])
    seqs_per_idx.update(dict(zip(subdf["orig_idx"], seqs_chr)))

# Reorder according to original index
seqs = [seqs_per_idx[i] for i in sorted(gr.df["orig_idx"])]

##
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/macaque/finetune"
model = keras.models.load_model(model_path, compile=False)

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

## Filter adata to specific peaks
# crested.pp.filter_regions_on_specificity(
#     adata, gini_std_threshold=1.0
# )

## ---------------------------
## Run enhancer sequences through model
crested.pl.heatmap.correlations_predictions(
    adata[consensus_anno["Group"].values, :],
    split="train",
    title="Correlations between Groundtruths and Predictions",
    x_label_rotation=90,
    width=15,
    height=15,
    log_transform=True,
    #reorder=False,
    vmax=1,
    vmin=-0.15,
    save_path=os.path.join(figure_dir, f"{species_name}_basalganglia_crested_correlations_predictions.png"),
)

crested.pl.heatmap.correlations_self(
    adata[consensus_anno["Group"].values, :],
    title="Self Correlation Heatmap",
    x_label_rotation=90,
    width=15,
    height=15,
    log_transform=True,
    #reorder=False,
    vmax=1,
    vmin=-0.15,
    save_path=os.path.join(figure_dir, f"{species_name}_basalganglia_crested_self_correlations.png"),
)

##
# plot predictions vs ground truth for a random region in the test set defined by index
test_df = adata.var[adata.var["split"] == "test"]
idx = 10000
region = test_df.index[idx]
print(region)
crested.pl.bar.region_predictions(adata[consensus_anno["Group"].values, :], 
                                    region, 
                                    title="Predictions vs Ground Truth", 
                                    x_label_rotation=90, 
                                    save_path=os.path.join(figure_dir, f"{species_name}_basalganglia_crested_region_{region}_predictions_vs_groundtruth.png"))