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
from pathlib import Path

out_path = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/human/ATAC/'
ocs_dir = '/allen/programs/celltypes/workgroups/rnaseqanalysis/bicore/OCS_multiome_datasets/human/'

## Get the task ID from the command line argument
task_id = int(sys.argv[1])-1

## Load in sample table
atac_samples = pd.read_csv(pjoin(out_path, 'atac_samples.csv'))
sample = atac_samples.iloc[task_id]

## Load in annotation table
anno_table = pd.read_csv(pjoin(out_path, 'anno_table.csv'), index_col=0)
anno_table['ar_cellname'] = anno_table['cell_barcode'].astype(str) + '-' + anno_table['ar_id'].astype(str)

## Process the sample
print(sample)
os.makedirs(out_path, exist_ok=True)
h5_out = pjoin(out_path, "h5ad", sample['ar_id'] + '.h5ad')

## Genome handling
genome_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/"
chr_sizes = pd.read_csv(pjoin(genome_path, 'star/chrNameLength.txt'), sep='\t', header=None)
chr_sizes_dict = dict(zip(chr_sizes[0], chr_sizes[1]))
fasta_path = pjoin(genome_path, 'fasta/genome.fa') if os.path.exists(pjoin(genome_path, 'fasta/genome.fa')) else pjoin(genome_path, 'fasta/genome.fa.gz')
gtf_path = pjoin(genome_path, 'genes/genes.gtf') if os.path.exists(pjoin(genome_path, 'genes/genes.gtf')) else pjoin(genome_path, 'genes/genes.gtf.gz')

genome = snap.genome.Genome(fasta=pjoin(genome_path, fasta_path),
                            annotation=pjoin(genome_path, gtf_path),
                            chrom_sizes=chr_sizes_dict)

## Load the fragment file
atac_h5ad = snap.pp.import_fragments(
    sample['fragment_file'],
    chrom_sizes=genome,
    file=h5_out,
    sorted_by_barcode=False,
)
atac_h5ad.obs['genome'] = [genome_path for x in range(atac_h5ad.shape[0])]
atac_h5ad.obs["cell_barcode"] = [re.sub('-[0-9]+','',x) for x in atac_h5ad.obs_names]
atac_h5ad.obs['ar_cellname'] = [re.sub('-[0-9]+','',x)+'-'+sample['ar_id'] for x in atac_h5ad.obs_names]

## Annotate
for anno in ['Neighborhood', 'Class', 'Subclass', 'Group','cluster_id', 'barcoded_cell_sample_label']:
    taxon_dict = dict(zip(anno_table['ar_cellname'], anno_table[anno]))
    atac_h5ad.obs[anno] = [str(taxon_dict.get(x,'nan')) for x in  atac_h5ad.obs['ar_cellname']]

## Subset to high quality RNA-seq cells
atac_h5ad.subset(obs_indices=np.array(atac_h5ad.obs_names)[atac_h5ad.obs['Group']!='nan'])

## Match cell ids with RNA-seq
atac_h5ad.obs_names = atac_h5ad.obs["cell_barcode"] + "-" + atac_h5ad.obs["barcoded_cell_sample_label"]

## Metrics / Features
snap.metrics.tsse(atac_h5ad, genome)
snap.pp.filter_cells(atac_h5ad, min_tsse=5.0, min_counts=1000)
snap.pp.add_tile_matrix(atac_h5ad, bin_size=5000)
snap.pp.select_features(atac_h5ad, n_features=50000)

## Clean up
atac_h5ad.close()

print(f"[DONE] Processed: {sample['ar_id']} → {h5_out}")
