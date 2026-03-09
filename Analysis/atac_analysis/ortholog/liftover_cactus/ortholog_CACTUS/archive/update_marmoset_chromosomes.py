import pandas as pd

##
peaks = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/marmoset/peaks/merged_peaks.bed", sep="\t", header=None)
peaks.columns = ["chrom", "start", "end"]

## alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

## Update chromosome names and leave if cannot map
peaks["chrom_revised"] = peaks["chrom"].map(dict(zip(chrom_alias["name"], chrom_alias["refseq"])))
peaks["chrom_revised"] = peaks["chrom_revised"].fillna(peaks["chrom"])
peaks["chrom_revised"].replace("v1_random", ".1", regex=True, inplace=True)

peaks["chrom"] = peaks["chrom_revised"]
peaks = peaks.drop(columns=["chrom_revised"])
peaks.to_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/marmoset/peaks/merged_peaks_updated_chrom.bed", sep="\t", header=False, index=False)