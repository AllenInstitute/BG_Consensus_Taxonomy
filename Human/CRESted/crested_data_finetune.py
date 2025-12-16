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
species_name = "human" # sys.argv[1] ## species_name = "mouse"
run_name = f"{species_name}_basalganglia_hmba_pre-print" ## Specifies a unique directory for this model

## Paths
work_dir = "./aibs-octo-dnaseq-modeling/basal-ganglia/"
model_dir = os.path.join(work_dir, "models", "crested")

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
adata = ad.read_h5ad(os.path.join(work_dir, "data", species_name, "crested_adata", f"{run_name}_crested.h5ad"))
adata_universe = adata.copy()

## ----------------------------------------------
## Specificity filtering 
crested.pp.filter_regions_on_specificity(
    adata, gini_std_threshold=1.0
)

## -----------------------------------------------
## Gene proximity filtering (avoiding promoters)

## Load the GTF file and extract transcript information
anno = pr.read_gtf(genome.annotation)
anno = anno[anno.Feature == "transcript"]

## Pick the transcript with the most upstream TSS per gene (by strand)
genes = anno.df.copy()
# genes.Chromosome = genes.Chromosome.replace(dict(zip(chromosomes["seqName"], chromosomes["UCSC"])))
genes["TSS"] = genes.apply(
    lambda row: row.Start if row.Strand == "+" else row.End, axis=1
)
genes_TSS = (
    genes
    .sort_values(["gene_id", "Strand", "TSS"])
    .groupby("gene_id", as_index=False)
    .first()
)

## Optional: Filter to a subset of gene names
genes_of_interest = ["DRD1", "DRD2", "STXBP6", "KCNIP1", "ARHGAP6", "PPP1R1B"]
genes_TSS["gene_id"] = genes_TSS["gene_name"].str.upper()
genes_TSS = genes_TSS[genes_TSS["gene_name"].isin(genes_of_interest)]

## Cis-regions: ±1 Mb
cis_df = genes_TSS.copy()
cis_df["Start"] = (cis_df["TSS"] - 1_000_000).clip(lower=0)
cis_df["End"] = cis_df["TSS"] + 1_000_000
cis_df = cis_df[["Chromosome", "Start", "End", "gene_name"]]
cis_ranges = pr.PyRanges(cis_df)

## Promoter regions: ±1 kb
promoter_df = genes_TSS.copy()
promoter_df["Start"] = (promoter_df["TSS"] - 1_000).clip(lower=0)
promoter_df["End"] = promoter_df["TSS"] + 1_000
promoter_df = promoter_df[["Chromosome", "Start", "End", "gene_name"]]
promoter_ranges = pr.PyRanges(promoter_df)

## Exclude promoter regions from cis-regions
cis_no_promoter = cis_ranges.subtract(promoter_ranges)

## Parse peak coordinates from adata.var.index = "chr:start-end"
peak_coords = adata_universe.var.index.to_series().str.extract(
    r'(?P<Chromosome>[^:]+):(?P<Start>\d+)-(?P<End>\d+)'
)
peak_coords.columns = ['Chromosome', 'Start', 'End']
peak_coords = peak_coords.astype({'Start': int, 'End': int})
peak_coords['region_name'] = adata_universe.var.index
peaks_pr = pr.PyRanges(peak_coords)

## Overlap peaks with cis regions excluding promoters
overlaps = peaks_pr.join(cis_no_promoter)
peaks_in_cis = overlaps.df['region_name'].unique()

## Fine tune set from specificity and proximity to marker genes
adata = adata_universe[:,adata_universe.var_names.isin(list(peaks_in_cis) + list(adata.var_names))].copy()

## Convert to long-form DataFrame
cell_types = adata.obs_names.tolist()
df = pd.DataFrame(adata.X.T, columns=cell_types)

##
medians = df.median()
q1 = df.quantile(0.25)
q3 = df.quantile(0.75)

plt.figure(figsize=(16, 6))
plt.bar(cell_types, medians, yerr=[medians - q1, q3 - medians])
plt.xticks(rotation=90)
plt.ylabel("Accessibility")
plt.title("Median Accessibility ± IQR per Cell Type")
plt.tight_layout()
plt.savefig(os.path.join("/home/nelson.johansen/", f"{run_name}_median_accessibility_per_cell_type.png"))

# --- Bar plot: how often each cell type is the most accessible per region ---
max_indices = adata.X.T.argmax(axis=1)  # index of max cell type per region
max_counts = np.bincount(max_indices, minlength=len(cell_types))

bar_df = pd.DataFrame({
    'CellType': cell_types,
    'MaxCount': max_counts
})

plt.figure(figsize=(16, 6))
sns.barplot(data=bar_df, x='CellType', y='MaxCount')
plt.xticks(rotation=90)
plt.title("Number of Regions Where Cell Type Has Max Accessibility")
plt.tight_layout()
plt.savefig(os.path.join("/home/nelson.johansen/", f"{run_name}_max_accessibility_per_cell_type.png"))

## Save finetuning adata object
adata.write_h5ad(os.path.join(work_dir, "data", species_name, "crested_adata", f"{run_name}_finetune_crested.h5ad"))