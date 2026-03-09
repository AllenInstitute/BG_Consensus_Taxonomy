import pandas as pd
import pyranges as pr
import os

## Directories
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## Species peak files
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
species_peak_files = {
    species: f"{data_dir}/{species}/ATAC/merged_peaks.bed"
    for species in ["human", "macaque", "marmoset"]
}

## Species genome annotation (GTF)
genome_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/ATAC/gencode"
species_gtf_files = {
    species: f"{genome_dir}/{species}_gencode.gtf.gz"
    for species in ["human", "macaque", "marmoset", "mouse"]
}

## Macaque chrom alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## Marmoset chrom alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

def parse_attributes(attributes_str):
    """Parse GTF attributes field into a dictionary."""
    attributes = {}
    for attr in attributes_str.strip().split(";"):
        if attr.strip():
            key_value = attr.strip().split(" ")
            if len(key_value) == 2:
                key = key_value[0]
                value = key_value[1].strip('"')
                attributes[key] = value
    return attributes

def build_promoter_regions(gtf_file, upstream_plus=1000, downstream_plus=500,
                           upstream_minus=500, downstream_minus=1000):
    """Return a PyRanges object of promoter regions from a GTF file."""
    gtf = pd.read_csv(gtf_file, sep="\t", header=None, comment="#",
                      names=["Chromosome", "Source", "Feature", "Start", "End",
                             "Score", "Strand", "Frame", "Attributes"], compression='gzip')
    if species == "macaque":
        ## Update for macaque chrom alias
        gtf["Chromosome"] = "chr" + gtf["Chromosome"].astype(str)
        gtf["Chromosome"].replace(dict(zip(macaque_chrom_alias[0], macaque_chrom_alias[2])), inplace=True)
    elif species == "marmoset":
        gtf["Chromosome"] = "chr" + gtf["Chromosome"].astype(str)
    ## Keep only gene features
    transcripts = gtf[gtf["Feature"] == "gene"].copy()
    ## Parse attributes
    transcripts["Attributes"] = transcripts["Attributes"].apply(parse_attributes)
    transcripts["gene_id"] = transcripts["Attributes"].apply(lambda x: x.get("gene_id", None))
    transcripts["gene_name"] = transcripts["Attributes"].apply(lambda x: x.get("gene_name", None))
    transcripts["gene_type"] = transcripts["Attributes"].apply(lambda x: x.get("gene_type", None))
    ## Keep protein-coding genes
    if transcripts["gene_type"].str.contains("protein_coding").any():
        transcripts = transcripts[transcripts["gene_type"] == "protein_coding"].copy()
    ## Strand-aware TSS
    tss = transcripts.copy()
    tss["TSS"] = tss.apply(lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1)  
    ## Expand to promoter region
    tss["Start"] = tss.apply(lambda r: max(r["TSS"] - upstream_plus, 0) if r["Strand"] == "+" else max(r["TSS"] - upstream_minus, 0), axis=1)
    tss["End"] = tss.apply(lambda r: r["TSS"] + downstream_plus if r["Strand"] == "+" else r["TSS"] + downstream_minus, axis=1)
    ## Create PyRanges object
    promoter_cols = ["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]
    return pr.PyRanges(tss[promoter_cols])

## Annotate peaks for all species
for species in ["marmoset"]: #species_peak_files.keys():
    print(f"Processing promoters for {species}...")
    ## Load peaks
    peaks = pr.read_bed(species_peak_files[species])
    ## Add Name to handle duplicate hits
    peaks = peaks.df
    peaks["Name"] = peaks["Chromosome"].astype(str) + ":" + peaks["Start"].astype(str) + "-" + peaks["End"].astype(str)
    peaks = pr.PyRanges(peaks)
    ## Build promoter regions
    promoters_gr = build_promoter_regions(species_gtf_files[species])
    ## Annotate peaks
    annotated = peaks.join(promoters_gr, how="left").df
    annotated = annotated.drop(columns=[c for c in annotated.columns if c.endswith("_b")])
    ## If gene_name is not null then assign promoter indicator
    annotated["promoter"] = annotated["gene_name"] != "-1"
    ## Keep the first
    annotated = annotated.drop_duplicates(subset=["Name"])
    if species == "marmoset":
        ## Map chrom aliases back to original
        annotated["Chromosome"] = annotated["Chromosome"].map(dict(zip(chrom_alias["name"], chrom_alias["refseq"])))
        annotated["Chromosome"] = annotated["Chromosome"].fillna(annotated["Name"].str.split(":").str[0])
    # Save results
    out_file = os.path.join(anno_dir, f"{species}_peaks_promoter.csv")
    annotated.to_csv(out_file, index=False)
    print(annotated.promoter.value_counts())
    print(f"Saved {out_file}")
