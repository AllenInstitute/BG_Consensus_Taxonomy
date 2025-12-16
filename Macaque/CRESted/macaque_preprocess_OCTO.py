import os
import sys
import crested
import keras
import glob
import re
import pandas as pd
import numpy as np

## -------------------------
## Parameters
species_name = "macaque"
run_name = "macaque_basalganglia_hmba"

## Paths
work_dir = os.path.abspath("./aibs-octo-dnaseq-modeling/")
model_dir = os.path.join(work_dir, "models", "crested")

## Genome handling
genome_path = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/{species_name}"
if species_name == "macaque":
    genome_path = os.path.join(genome_path, "ncbi")

##
chr_sizes = os.path.join(genome_path, 'star/chrNameLength.txt')
fasta_path = os.path.join(genome_path, 'fasta/genome.fa') if os.path.exists(os.path.join(genome_path, 'fasta/genome.fa')) else os.path.join(genome_path, 'fasta/genome.fa.gz')
gtf_path = os.path.join(genome_path, 'genes/genes.gtf') if os.path.exists(os.path.join(genome_path, 'genes/genes.gtf')) else os.path.join(genome_path, 'genes/genes.gtf.gz')

## -------------------------
bigwigs_folder = os.path.join(work_dir, "basal-ganglia", "data", species_name, "bigwigs")
regions_file = os.path.join(work_dir, "basal-ganglia", "data", species_name, "peaks", f"merged_peaks.bed")

## Load in the bigwigs and build the crested anndata object
adata = crested.import_bigwigs(
    bigwigs_folder=bigwigs_folder,
    regions_file=regions_file,
    chromsizes_file=chr_sizes,
    target_region_width=1000, ## Standard peak window from MACS3 (peak-caller)
    target="count",
)
adata.X = np.nan_to_num(adata.X, nan=0.0)

## Extend Seqeunce to 2114 but only predicting the center 500bp ATAC-seq profile as defined above.
crested.pp.change_regions_width(
    adata,
    width=2114,
    chromsizes_file=chr_sizes,
)

##
crested.pp.normalize_peaks(
    adata, top_k_percent=0.03
)

crested.pl.bar.normalization_weights(
    adata, title="Normalization Weights per Cell Type", x_label_rotation=90, height=12, width=16, save_path=f"/home/nelson.johansen/TSS_norm_weights_{species_name}.png"
)

## Train, validation, test split
crested.pp.train_val_test_split(
    adata, strategy="chr", val_chroms=["NC_041763.1"], test_chroms=["NC_041771.1"]
)

## ------------------------------------------------------
## Filter to only full chromosomes and remove contigs
chrom_alias = pd.read_table(os.path.join(genome_path, "rheMac10.chromAlias.txt"), delimiter="\t")
chrom_alias.columns = ["seqName", "aliasName", "UCSC"]

## Gather chr1-20, chrX,Y
chromosomes = chrom_alias.loc[chrom_alias.UCSC.str.startswith("NC_")]
chromosomes = chromosomes.loc[~chromosomes.seqName.str.contains("chrM")]

## Keep only desired chromosomes and filter away contigs
adata = adata[:,adata.var.index[adata.var.chr.isin(chromosomes.UCSC)]].copy()

## You will want to save this file for evaluation / inference steps.
if ~os.path.exists(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata")):
    os.makedirs(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata"), exist_ok=True)

## Write filtered file
adata.write_h5ad(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata", f"{run_name}_pre-print_crested.h5ad"))
