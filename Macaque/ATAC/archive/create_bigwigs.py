import snapatac2 as snap
import anndata as ad
import pandas as pd
import glob
import os
import sys

## Where mapped data lives
mapped_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Macaque/"

## Where your ATAC concatenated h5ads are
os.chdir("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC/h5ads")

try:
    print("loading subset object")
    adata_subset = snap.read_dataset("macaque_BG_AIT117_subset/_dataset.h5ads")
except:
    print("rebuilding subset")
    ## Load in ATAC
    adata = snap.read_dataset("macaque_BG_AIT117.h5ads")
    adata.obs["batch"] = [batch.split("/")[-1].replace(".h5ad", "") for batch in adata.obs["sample"].to_list()]
    adata.obs_names = [f"{obs_name}-{batch}" for obs_name, batch in zip(adata.obs_names, adata.obs["batch"])]
    adata.obs_names = ["-".join(obs_name.split("-")[0:3]) for obs_name in adata.obs_names] ## This is here in case you load the file multiple times. 
    ## Load in mapped data
    mapping_file = os.path.join(mapped_dir, "HMBA_Macaque_BG_082024_AIT.h5ad")
    mapping = ad.read_h5ad(mapping_file, backed="r")
    mapping.obs["barcode"] = [barcode.split("-")[0] for barcode in mapping.obs.cell_id]
    mapping.obs["Sample"] = [sample.split("-")[-1] for sample in mapping.obs.index]
    mapping.obs.index = mapping.obs["barcode"] + "-1-macaque_" + mapping.obs.Sample.astype(str) 
    ## Find the barcodes / cell_ids in common between RNA and ATAC
    common_cells = list(set(adata.obs_names) & set(mapping.obs.index))
    ## Build a subset of the concatenated ATAC data, this creates a new file and directory
    adata_subset = adata.subset(common_cells, out="macaque_BG_AIT117.h5ads".replace(".h5ads", "_subset"))[0]
    adata.close() ## Close the full ATAC file as we don't need it anymore
    ## Neighborhood, Class, Subclass, Group annotations
    adata_subset.obs["Neighborhood"] = mapping.obs.loc[adata_subset.obs_names, "Neighborhood_label"].str.replace(" ", "_")
    adata_subset.obs["Class"] = mapping.obs.loc[adata_subset.obs_names, "Class_label"].str.replace(" ", "_")
    adata_subset.obs["Subclass"] = mapping.obs.loc[adata_subset.obs_names, "Subclass_label"].str.replace(" ", "_")
    adata_subset.obs["Group"] = mapping.obs.loc[adata_subset.obs_names, "Group_label"].str.replace(" ", "_")

##
# bigwig_dir = "/allen/programs/celltypes/workgroups/hct/NelsonJ/group"
# os.makedirs(bigwig_dir, exist_ok=True)

# ##
# snap.ex.export_coverage(adata_subset, groupby='Group', normalization="CPM", counting_strategy="insertion", suffix='.bw', out_dir=bigwig_dir, n_jobs=4)

##################################################################
## Call peaks and build consensus set across annotations
snap.tl.macs3(adata_subset, groupby='Group', n_jobs=24)

## -------
annot_tab=pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Macaque/AIT117/ATAC/ARC_refs_tab.csv", index_col=0)
references = dict(zip(annot_tab['species'], annot_tab['annotation']))
reference = references["macaque"]

## Chromosome sizes
chr_sizes = pd.read_csv(os.path.join(reference, 
                                        'star/chrNameLength.txt'), sep='\t', header=None)
chr_sizes_dict = dict(zip(chr_sizes[0], chr_sizes[1]))

## Fasta and GTF paths
fasta_path = os.path.join(reference, 'fasta/genome.fa') if os.path.exists(os.path.join(reference, 'fasta/genome.fa')) else os.path.join(reference, 'fasta/genome.fa.gz')
gtf_path = os.path.join(reference, 'genes/genes.gtf') if os.path.exists(os.path.join(reference, 'genes/genes.gtf')) else os.path.join(reference, 'genes/genes.gtf.gz')

## -------------------------
genome = snap.genome.Genome(fasta=os.path.join(reference, fasta_path),
                            annotation=os.path.join(reference, gtf_path),
                            chrom_sizes=chr_sizes_dict)

## Call consensus peaks
consensus_peaks = snap.tl.merge_peaks(adata_subset.uns['macs3'], chrom_sizes=genome)

## Setup bed file
consensus_peaks_pd = consensus_peaks.to_pandas()
consensus_peaks_pd["chr"] = [peak.split(":")[0] for peak in consensus_peaks_pd["Peaks"]]
consensus_peaks_pd["start"] = [peak.split(":")[1].split("-")[0] for peak in consensus_peaks_pd["Peaks"]]
consensus_peaks_pd["end"] = [peak.split(":")[1].split("-")[1] for peak in consensus_peaks_pd["Peaks"]]

# Write to a BED file
peak_dir = os.path.join("/allen/programs/celltypes/workgroups/hct/NelsonJ/group")
os.makedirs(peak_dir, exist_ok=True)

with open(os.path.join(peak_dir, "hmba_macaque_AIT117_consensus_peaks.bed"), "w") as f:
    consensus_peaks_pd[["chr", "start", "end"]].to_csv(f, sep="\t", header=False, index=False)

with open(os.path.join(peak_dir, "hmba_macaque_AIT117_consensus_peaks_celltypes.bed"), "w") as f:
    consensus_peaks_pd.to_csv(f, sep="\t", header=True, index=False)

##
adata_subset.close()