import os, sys
import anndata
import pandas as pd
import snapatac2 as snap

## -------------------------
## Get the task ID from the command line argument, SLURM task array
task_id = int(sys.argv[1])
print(task_id)

## -------------------------
## Where are we going to save all the output files
wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC"
os.chdir(wkdir)

## -------------------------
## Load the reference annotation table, essentially a look up for genome files (fasta, gtf) handled by the BICore
annot_tab=pd.read_csv("ARC_refs_tab.csv", index_col=0)
references = dict(zip(annot_tab['species'], annot_tab['annotation']))

## Load in LIMS tracker file for M1 ATAC
ATAC_tracker = pd.read_csv("ATAC_HMBA_Macaque_BG_study.csv")
ATAC_tracker.loc[ATAC_tracker['species'] == 'Macaca nemestrina', 'species'] = 'macaque'
ATAC_tracker.loc[ATAC_tracker['species'] == 'Macaca mulatta', 'species'] = 'macaque'

##
ATAC_tracker['reference'] = ATAC_tracker['species'].replace(references)

## Create the ATAC samples DataFrame
ATAC_tracker['path'] = ATAC_tracker.ar_directory + ATAC_tracker.ar_id
ATAC_tracker['fragment_file'] = ATAC_tracker['path'] + '/outs/atac_fragments.tsv.gz'

## Gather all information for atac samples
atac_samples = ATAC_tracker.loc[:, ['species', 'fragment_file', 'reference', 'path', 'ar_directory', 'ar_id']]
sample = atac_samples.iloc[task_id]

## -------------------------
## Process the sample
print(sample)
os.makedirs(os.path.join(wkdir, sample.ar_id), exist_ok=True)

## Chromosome sizes
chr_sizes = pd.read_csv(os.path.join(sample['reference'], 
                                        'star/chrNameLength.txt'), sep='\t', header=None)
chr_sizes_dict = dict(zip(chr_sizes[0], chr_sizes[1]))

## Fasta and GTF paths
fasta_path = os.path.join(sample['reference'], 'fasta/genome.fa') if os.path.exists(os.path.join(sample['reference'], 'fasta/genome.fa')) else os.path.join(sample['reference'], 'fasta/genome.fa.gz')
gtf_path = os.path.join(sample['reference'], 'genes/genes.gtf') if os.path.exists(os.path.join(sample['reference'], 'genes/genes.gtf')) else os.path.join(sample['reference'], 'genes/genes.gtf.gz')

## -------------------------
genome = snap.genome.Genome(fasta=os.path.join(sample['reference'], fasta_path),
                            annotation=os.path.join(sample['reference'], gtf_path),
                            chrom_sizes=chr_sizes_dict)

## Output file
h5_out = os.path.join(wkdir, sample.ar_id, sample['species'] + "_" + sample['ar_id'] + '.h5ad')

## -------------------------
# if os.path.exists(h5_out):
#     print(f"Output file {h5_out} already exists. Skipping.")
#     sys.exit(0)

## -------------------------
## Import the fragment file and process using snapatac2
print(genome)
print(sample['reference'])
adata_atac = snap.pp.import_data(
    sample['fragment_file'],
    chrom_sizes=genome,
    file=h5_out,
    sorted_by_barcode = False,
)

## Add genome
adata_atac.uns['genome'] = sample['reference']
adata_atac.uns['species'] = sample['species']
adata_atac.uns['ar_id'] = sample['ar_id']

## -------------------------
## Add metrics 
snap.metrics.tsse(adata_atac, genome)
snap.metrics.frag_size_distr(adata_atac)

## Close file
adata_atac.close()