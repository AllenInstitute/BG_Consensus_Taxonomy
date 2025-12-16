import snapatac2 as snap
import anndata as ad
import pandas as pd
import glob
import os
import sys

##
os.chdir("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC/h5ads")

## Load in subset data with labels
adata = snap.read_dataset("macaque_subset/_dataset.h5ads")

## Call peaks and build consensus set across annotations
snap.tl.macs3(adata, groupby='Group')

## -------
annot_tab=pd.read_csv("ARC_refs_tab.csv", index_col=0)
references = dict(zip(annot_tab['species'], annot_tab['annotation']))
reference = references["macaque"]

## Chromosome sizes
chr_sizes = pd.read_csv(os.path.join(reference, 
                                        'mm10.chrom.sizes'), sep='\t', header=None)
chr_sizes_dict = dict(zip(chr_sizes[0], chr_sizes[1]))

## Fasta and GTF paths
fasta_path = os.path.join(reference, 'fasta/mm10.fa')
gtf_path = os.path.join(reference, 'genes/mm10.ncbiRefSeq.gtf.gz')

## -------------------------
genome = snap.genome.Genome(fasta=os.path.join(reference, fasta_path),
                            annotation=os.path.join(reference, gtf_path),
                            chrom_sizes=chr_sizes_dict)

## Call consensus peaks
consensus_peaks = snap.tl.merge_peaks(adata.uns['macs3'], chrom_sizes=genome)

adata.close()

## Setup bed file
consensus_peaks_pd = consensus_peaks.to_pandas()
consensus_peaks_pd["chr"] = [peak.split(":")[0] for peak in consensus_peaks_pd["Peaks"]]
consensus_peaks_pd["start"] = [peak.split(":")[1].split("-")[0] for peak in consensus_peaks_pd["Peaks"]]
consensus_peaks_pd["end"] = [peak.split(":")[1].split("-")[1] for peak in consensus_peaks_pd["Peaks"]]

# Write to a BED file
peak_dir = os.path.join("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC/peaks", species_name)
os.makedirs(peak_dir, exist_ok=True)

with open(os.path.join(peak_dir, "hmba_macaque_AIT117_consensus_peaks.bed"), "w") as f:
    consensus_peaks_pd[["chr", "start", "end"]].to_csv(f, sep="\t", header=False, index=False)
