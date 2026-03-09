import pandas as pd
from tqdm  import tqdm
import pyranges as pr
import os

## Directories
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## Column names based on UCSC schema for rmsk
colnames = [
    "bin",        # UCSC binning index (for fast overlap queries)
    "swScore",    # Smith-Waterman alignment score
    "milliDiv",   # divergence from consensus (per 1000 bases)
    "milliDel",   # deletions relative to consensus (per 1000 bases)
    "milliIns",   # insertions relative to consensus (per 1000 bases)
    "Chromosome",   # chromosome/contig name
    "Start",        # alignment start on the genome (0-based)
    "End",          # alignment end on the genome
    "genoLeft",   # remaining bases in chromosome after alignment
    "Strand",     # orientation of the repeat (+ or -)
    "repName",    # name of the repeat (e.g., AluY, L1HS)
    "repClass",   # repeat class (e.g., SINE, LINE, LTR, DNA, Simple_repeat)
    "repFamily",  # repeat family (e.g., Alu, L1, hAT-Charlie)
    "repStart",   # start in the repeat consensus sequence
    "repEnd",     # end in the repeat consensus sequence
    "repLeft",    # remaining bases in the repeat consensus
    "id"          # unique identifier
]

## Peak data
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
species_peak_files = {
    #"human": os.path.join(data_dir, "human/ATAC/merged_peaks.bed"),
    #"macaque": os.path.join(data_dir, "macaque/ATAC/merged_peaks.bed"),
    "marmoset": os.path.join(data_dir, "marmoset/ATAC/merged_peaks.bed")
}

## Identify repeat masker annotation file
genome_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/ATAC/"
species_rmsk_files = {
    species: os.path.join(genome_dir, f"rmsk_{species}_goldenUCSC.txt.gz")
    for species in ["hg38", "rheMac10"]
}
species_rmsk_files["calJacX"] = os.path.join(genome_dir, "GCF_011100555.1.repeatMasker.out.gz")

## Species name to genome assembly
species_to_assembly = {"human": "hg38", "macaque": "rheMac10", "marmoset": "calJacX"}

## Macaque chrom alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## Marmoset chrom alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

## Annotate peaks with TEs
for species in species_peak_files.keys():
    print(f"Processing {species}...")
    ## Read in peaks
    species_bed = pr.read_bed(species_peak_files[species])
    ## Add Name to handle duplicate hits
    species_bed = species_bed.df
    species_bed["Name"] = species_bed["Chromosome"].astype(str) + ":" + species_bed["Start"].astype(str) + "-" + species_bed["End"].astype(str)
    species_bed = pr.PyRanges(species_bed)
    ## Read in rmsk
    if species == "macaque":
        rmsk = pd.read_csv(species_rmsk_files[species_to_assembly[species]], sep="\t", header=None, names=colnames, compression='gzip')
        ## Update chromosome names
        rmsk["Chromosome"] = rmsk["Chromosome"].astype(str)
        rmsk["Chromosome"].replace(dict(zip(macaque_chrom_alias[0], macaque_chrom_alias[2])), inplace=True)
    elif species == "marmoset":
        rmsk = pd.read_csv(species_rmsk_files[species_to_assembly[species]], sep="\s+", compression='gzip', skiprows=1, index_col=False)
        rmsk.rename(columns={
            "score": "swScore",
            "div.": "milliDiv",
            "del.": "milliDel",
            "ins.": "milliIns",
            "sequence": "Chromosome",
            "begin": "Start",
            "end": "End",
            "(left)": "genoLeft",
            "class/family": "repName",
            "ID": "id"
        }, inplace=True)
        rmsk["repClass"] = rmsk["begin.1"].apply(lambda x: x.split("/")[0])
        rmsk["repFamily"] = rmsk["begin.1"].apply(lambda x: x.split("/")[1] if "/" in x else "-1")
        ## Update chromosome names
        rmsk["Chromosome"] = rmsk["Chromosome"].map(dict(zip(chrom_alias["refseq"], chrom_alias["name"])))
    else:
        rmsk = pd.read_csv(species_rmsk_files[species_to_assembly[species]], sep="\t", header=None, names=colnames, compression='gzip')
    ## Convert to PyRanges
    TEs = pr.PyRanges(rmsk)
    ## Annotate peaks with TEs
    annotated = species_bed.join(TEs, how="left", report_overlap=True).df
    annotated = annotated.drop(columns=[c for c in annotated.columns if c.endswith("_b")])
    ## Remove results for overlap < 250bp (50% of peak)
    annotated.loc[annotated["Overlap"] < 250, ["repName", "repClass", "repFamily"]] = "-1"
    ## If repClass is of LINE, SINE, LTR, DNA or RC then assign TE indicator
    annotated["TE"] = annotated["repClass"].apply(lambda x: True if x in ["LINE", "SINE", "LTR", "DNA", "RC"] else False)
    ## Summarize overlapping TEs per peak
    summary = annotated.groupby("Name").agg({
        "TE": "any",  ## True if any overlap
        "repName": lambda x: ",".join(x.dropna().unique()),
        "repClass": lambda x: ",".join(x.dropna().unique()),
        "repFamily": lambda x: ",".join(x.dropna().unique())
    }).reset_index()
    ## Merge back to annotated
    annotated_TE = species_bed.df.merge(summary, on="Name", how="left")
    annotated_TE.fillna({"TE": False, "repName": "-1", "repClass": "-1", "repFamily": "-1"}, inplace=True)
    ## Map chrom aliases back to original for marmoset
    if species == "marmoset":
        ## Map chrom aliases back to original
        annotated_TE["Chromosome"] = annotated_TE["Chromosome"].map(dict(zip(chrom_alias["name"], chrom_alias["refseq"])))
        annotated_TE["Chromosome"] = annotated_TE["Chromosome"].fillna(annotated_TE["Name"].str.split(":").str[0])
    print(annotated_TE.head())
    ## Print stats
    print(f"%TE: {annotated_TE['TE'].mean()*100:.2f}% ({annotated_TE['TE'].sum()} / {annotated_TE.shape[0]})")
    ## Save results
    annotated_TE.to_csv(os.path.join(anno_dir, f"{species}_peaks_annotated_TE.csv"), index=False)

