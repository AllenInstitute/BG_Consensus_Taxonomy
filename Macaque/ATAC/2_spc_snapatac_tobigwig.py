import os
import pandas as pd
import re
import sys
import matplotlib.pyplot as plt
import scipy
import numpy as np
import anndata
import h5py
import tqdm
from os.path import join as pjoin
import snapatac2 as snap

root_directory = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/ATAC'
groupby = 'Group'

##
file_path = os.path.join(root_directory, 'concatenated.h5ads')
data = snap.read_dataset(file_path, mode='r+')
data.obs[groupby] = data.adatas.obs[groupby]

## bigwigs directorys
out_dir=os.path.join(os.path.dirname(root_directory), f'{groupby}_bigwig_TSSnorm')
os.makedirs(out_dir, exist_ok=True)

## Genome handling
genome_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi"
chr_sizes = pd.read_csv(pjoin(genome_path, 'star/chrNameLength.txt'), sep='\t', header=None)
chr_sizes_dict = dict(zip(chr_sizes[0], chr_sizes[1]))
fasta_path = pjoin(genome_path, 'fasta/genome.fa') if os.path.exists(pjoin(genome_path, 'fasta/genome.fa')) else pjoin(genome_path, 'fasta/genome.fa.gz')
gtf_path = pjoin(genome_path, 'genes/genes.gtf') if os.path.exists(pjoin(genome_path, 'genes/genes.gtf')) else pjoin(genome_path, 'genes/genes.gtf.gz')

genome = snap.genome.Genome(fasta=pjoin(genome_path, fasta_path),
                            annotation=pjoin(genome_path, gtf_path),
                            chrom_sizes=chr_sizes_dict)

gtf = pybedtools.BedTool(genome.annotation)

# Filter for 'transcript' features
transcripts = gtf.filter(lambda x: x[2] == "transcript")

# Function to generate ±100bp window around TSS
def tss_window(feature):
    if feature.strand == "+":
        tss = feature.start  # 0-based for BED
        start = max(tss - 100, 0)  # ensure non-negative start
        end = tss + 101  # end is exclusive
    else:
        tss = feature.end  # 0-based BED end is exclusive
        start = max(tss - 101, 0)
        end = tss + 100
    gene_id = feature.attrs.get("gene_id", "NA")
    return (feature.chrom, start, end, gene_id, 0, feature.strand)

# Create list of TSS ±100bp intervals
tss_windows = [tss_window(f) for f in transcripts]

# Convert to DataFrame
tss_df = pd.DataFrame(tss_windows, columns=["chrom", "start", "end", "gene_id", "score", "strand"])
tss_df = tss_df.drop_duplicates()
# Convert each row to "chrom:start-end" format
tss_coords = tss_df.apply(lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1)
# Convert to list of strings
tss_list = tss_coords.tolist()

##
snap.ex.export_coverage(data, groupby=groupby,out_dir=out_dir, normalization='CPM', bin_size=10, counting_strategy='insertion', include_for_norm=tss_list)

# ## Compute bigwig
# snap.ex.export_coverage(data, groupby=groupby,out_dir=out_dir,normalization='CPM', bin_size=10, counting_strategy='insertion')

# ## Call peaks
# snap.tl.macs3(data,groupby=groupby) 

# ## Export peaks
# peak_tables={}
# for k in data.uns['macs3'].keys():
#     if 'dict' in str(type(data.uns['macs3'][k])):
#         peak_tables[k]=list(data.uns['macs3'][k].values())[0]
#     else:
#         peak_tables[k]=data.uns['macs3'][k]

# df = snap.tl.merge_peaks(peak_tables,
#                         dict(zip(data.uns['reference_sequences']['reference_seq_name'],data.uns['reference_sequences']['reference_seq_length'])))
# df.write_csv(os.path.join(os.path.dirname(file_path),f'{groupby}_bigwig',f"{groupby}_by_peaks.csv"))
# merged_df=pd.DataFrame(list(df['Peaks'].to_pandas().str.split('\:|\-')),columns=['chrom','start','end'])
# merged_df.name='.'
# merged_df.to_csv(os.path.join(os.path.dirname(file_path),f'{groupby}_bigwig',"merged_peaks.bed"),sep='\t',header=None,index=False)
data.close()