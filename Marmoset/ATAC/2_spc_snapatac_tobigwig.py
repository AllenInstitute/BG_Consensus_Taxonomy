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
import anndata as ad
import pybedtools

root_directory = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/marmoset/ATAC'
os.chdir(root_directory)
groupby = 'Group'

try:
    print("loading subset object")
    adata_subset = snap.read_dataset("_subset/_dataset.h5ads")
except:
    ##
    file_path = os.path.join(root_directory, 'concatenated.h5ads')
    data = snap.read_dataset(file_path, mode='r+')
    ## 
    # alignment_id = [f.split("_")[0] for f in data.adatas.obs["sample_name"]]
    # alignment_id_series = pd.Series(alignment_id, index=data.obs_names)
    # obs_names_series = pd.Series(data.obs_names, index=data.obs_names)
    # data.obs["sample_id"] = obs_names_series + "_" + alignment_id_series
    ##
    rna_directory = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/s3/AIT'
    file_path = os.path.join(rna_directory, 'Marmoset_HMBA_basalganglia_AIT_pre-print.h5ad')
    rna_taxonomy = ad.read_h5ad(file_path, backed='r')
    # rna_taxonomy.obs["sample_id"] = rna_taxonomy.obs.cell_barcode.astype(str) + "-1" + "_" + rna_taxonomy.obs.alignment_job_id.astype(str)
    ## Merge rna and atac
    common_cells = list(set(data.obs_names).intersection(set(rna_taxonomy.obs_names)))
    ## Build a subset of the concatenated ATAC data, this creates a new file and directory
    adata_subset = data.subset(common_cells, out="_subset")[0]
    data.close() ## Close the full ATAC file as we don't need it anymore
    ## Neighborhood, Class, Subclass, Group annotations
    adata_subset.obs["Neighborhood"] = rna_taxonomy.obs.loc[adata_subset.obs_names, "Neighborhood"].str.replace(" ", "_")
    adata_subset.obs["Class"] = rna_taxonomy.obs.loc[adata_subset.obs_names, "Class"].str.replace(" ", "_")
    adata_subset.obs["Subclass"] = rna_taxonomy.obs.loc[adata_subset.obs_names, "Subclass"].str.replace(" ", "_")
    adata_subset.obs["Group"] = rna_taxonomy.obs.loc[adata_subset.obs_names, "Group"].str.replace(" ", "_")

## bigwigs directorys
out_dir=os.path.join(os.path.dirname(root_directory), f'{groupby}_bigwig_TSSnorm')
os.makedirs(out_dir, exist_ok=True)

## Genome handling
genome_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/marmoset/hmba"
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
snap.ex.export_coverage(adata_subset, groupby=groupby,out_dir=out_dir, normalization='CPM', bin_size=10, counting_strategy='insertion', include_for_norm=tss_list, n_jobs=24)

## Call peaks
snap.tl.macs3(adata_subset, groupby=groupby, n_jobs=24) 

## Export peaks
peak_tables={}
for k in adata_subset.uns['macs3'].keys():
    if 'dict' in str(type(adata_subset.uns['macs3'][k])):
        peak_tables[k]=list(adata_subset.uns['macs3'][k].values())[0]
    else:
        peak_tables[k]=adata_subset.uns['macs3'][k]

df = snap.tl.merge_peaks(peak_tables,
                        dict(zip(adata_subset.uns['reference_sequences']['reference_seq_name'],adata_subset.uns['reference_sequences']['reference_seq_length'])))
df.write_csv(os.path.join(os.path.dirname(out_dir),f"{groupby}_by_peaks.bed"))
merged_df=pd.DataFrame(list(df['Peaks'].to_pandas().str.split('\:|\-')),columns=['chrom','start','end'])
merged_df.name='.'
merged_df.to_csv(os.path.join(os.path.dirname(out_dir),"merged_peaks.bed"),sep='\t',header=None,index=False)
adata_subset.close()