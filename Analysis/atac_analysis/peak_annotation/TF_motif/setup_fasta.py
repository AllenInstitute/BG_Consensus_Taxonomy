import os, sys
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
from tqdm import tqdm
import pyranges as pr
from pyfaidx import Fasta

## Function to fetch sequences in batch
def fetch_sequences(chrom, starts, ends):
    ref = fasta[chrom]
    return [ref[start:end].seq for start, end in zip(starts, ends)]

## Function to fetch sequences in order
def fetch_sequences(gr_df, fasta, fetch_fn, seq_len=500):
    """
    Fetch sequences from a FASTA file for regions in a PyRanges dataframe,
    preserving the original order.

    Parameters
    ----------
    gr_df : pandas.DataFrame
        A dataframe with columns ['Chromosome', 'Start', 'End'] at minimum.
    fasta : Fasta
        An open pyfaidx.Fasta object.
    fetch_fn : function
        A function with signature (chrom, starts, ends) -> list of sequences.
    seq_len : int, optional
        Length of sequence to insert if chromosome not found. Default 500.

    Returns
    -------
    list of str
        Sequences in the same order as the original gr_df.
    """
    gr_df = gr.copy()
    gr_df["orig_idx"] = gr_df.index
    seqs_per_idx = {}
    for chrom, subdf in tqdm(gr_df.groupby("Chromosome"), desc="Fetching sequences"):
        if chrom not in fasta.keys():
            seqs_chr = ["N" * seq_len] * len(subdf)
        else:
            seqs_chr = fetch_fn(chrom, subdf["Start"], subdf["End"])
        seqs_per_idx.update(dict(zip(subdf["orig_idx"], seqs_chr)))
    ## Reorder to match original index
    ordered_seqs = [seqs_per_idx[i] for i in sorted(gr_df["orig_idx"])]
    return ordered_seqs

## Load bed file
species = "human"
bed_file = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/{species}/peaks/merged_peaks.bed"

bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "end"])

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

## Convert BED to PyRanges
gr = pr.PyRanges(bed_df.rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'}))

## Function to fetch sequences in batch
seqs = fetch_sequences(gr, fasta, fetch_sequences, seq_len=500)

## Write out FASTA file give each sequence a chr:start-end header
output_fasta = os.path.dirname(bed_file) + "/merged_peaks.fasta"
with open(output_fasta, "w") as f:
    for i, row in tqdm(bed_df.iterrows()):
        header = f">{row['chr']}:{row['start']}-{row['end']}\n"
        sequence = seqs[i] + "\n"
        f.write(header)
        f.write(sequence)


## grep -h '^MOTIF' -A 10000 path/to/meme_files/*.meme >> merged.meme