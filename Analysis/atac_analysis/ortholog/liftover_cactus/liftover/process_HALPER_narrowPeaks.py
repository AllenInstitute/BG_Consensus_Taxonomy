import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import gzip
import glob 
import os

import pyranges as pr
import pandas as pd

## Peak data
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
species_peaks = {
    "Homo_sapiens": pd.read_csv(os.path.join(data_dir, "human/ATAC/merged_peaks.bed"), 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
    "Macaca_mulatta": pd.read_csv(os.path.join(data_dir, "macaque/ATAC/merged_peaks.bed"), 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
    "Callithrix_jacchus": pd.read_csv(os.path.join(data_dir, "marmoset/ATAC/merged_peaks.bed"), 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None),
}

## halLiftover data
liftover_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/"

## Narrowpeak file columns and read
names = ["Chromosome", "Start", "End", "summit_pos", "ID", "target_length", "query_length", "target_summit_to_target_ortholog_start_length", "target_summit_to_target_ortholog_end_length"]
species_liftovers = {
    "Macaca_mulatta": pd.read_csv(os.path.join(liftover_dir, "macaque/merged_peaks.Macaca_mulattaToHomo_sapiens.HALPER_.narrowPeak.gz"), 
                                    names=names, sep="\t", header=None, compression="gzip"),
    "Callithrix_jacchus": pd.read_csv(os.path.join(liftover_dir, "marmoset/merged_peaks_updated_chrom.Callithrix_jacchusToHomo_sapiens.HALPER_.narrowPeak.gz"), 
                                  names=names, sep="\t", header=None, compression="gzip")
}

## alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

## Now merge into liftover set anchoring on human when possible
def build_liftover_set(human_bed, species_liftovers, minMatch=250):
    """
    Build human-anchored liftover peak set.

    Parameters
    ----------
    human_bed : pd.DataFrame
        Human peaks with columns ["Chromosome","Start","End","ID"].
    species_liftovers : dict[str, pd.DataFrame]
        Dict mapping species -> liftover peaks in human coords
        (each with ["Chromosome","Start","End","ID"]).

    Returns
    -------
    pd.DataFrame
        Ortholog peak set with columns:
        ["ortholog_region_id","Chromosome","Start","End","source","ID"]
    """
    ## Anchor on human peaks
    human_df = human_bed.copy()
    human_df["ID"] = human_df["Chromosome"].astype(str) + ":" + human_df["Start"].astype(str) + "-" + human_df["End"].astype(str)
    ## Set index without giving a name
    human_df.set_index("ID", inplace=True, drop=False)
    ## Convert human peaks to PyRanges
    human_pr = pr.PyRanges(human_df)
    for species, df in species_liftovers.items():
        print(f"Processing {species}...")
        ## Remove rows without chromosome and with small overlap
        df = df[df["Chromosome"] != ""].copy()
        df = df[df.target_length >= minMatch].copy() ## Only attempt merge in if >=50% of bases aligned to Human
        if species == "Callithrix_jacchus":
            df["ID"] = df["ID"].str.replace(r':\d+$', '', regex=True) ## Correction for allowing HALPER to set peak ids
            df["Chromosome-species"] = df["ID"].str.split(":").str[0]
            df["Start-species"] = df["ID"].str.split(":").str[1].str.split("-").str[0].astype(int)
            df["End-species"] = df["ID"].str.split(":").str[1].str.split("-").str[1].astype(int)
            df["Chromosome-species"] = df["Chromosome-species"].map(dict(zip(chrom_alias["refseq"], chrom_alias["name"])))
            df["ID"] = df["Chromosome-species"].astype(str) + ":" + df["Start-species"].astype(str) + "-" + df["End-species"].astype(str)
        if species == "Macaca_mulatta":
            df["ID"] = df["ID"].str.replace(r':\d+$', '', regex=True) ## Correction for allowing HALPER to set peak ids
        ## Assign alignment ID
        df["alignment_id"] = df["Chromosome"] + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
        ## Convert to pyrange
        species_pr = pr.PyRanges(df.rename(columns={"ID":f"{species}_ID"}))
        ## Find overlaps with human peaks
        intersect_df = human_pr.join(species_pr).as_df()
        #intersect_df = human_pr.join(species_pr, how="right").as_df()
        #intersect_df = intersect_df[intersect_df["ID"] != "-1"].copy()
        human_df.loc[intersect_df.ID.tolist(), f"{species}_ID"] = intersect_df[f"{species}_ID"].values
        human_df.loc[intersect_df.ID.tolist(), f"{species}_aligned_ID"] = intersect_df["alignment_id"].values
    ## Rename human ID column
    human_df.rename(columns={"ID":"human_ID"}, inplace=True)
    ## Return
    return human_df

## Create human-anchored liftover set
liftover_set = build_liftover_set(
    species_peaks["Homo_sapiens"],
    species_liftovers,
    minMatch=250 ## 50% overlap for 500bp peaks
)

## Identify human specific peaks by is.na on all but human_ID columns
liftover_set["human_specific"] = False
liftover_set.loc[
    (liftover_set["Macaca_mulatta_ID"].isna()) & 
    (liftover_set["Callithrix_jacchus_ID"].isna()),
    "human_specific"
] = True

liftover_set["ortholog"] = False
liftover_set.loc[
    (~liftover_set["Macaca_mulatta_ID"].isna()) & 
    (~liftover_set["Callithrix_jacchus_ID"].isna()),
    "ortholog"
] = True

liftover_set.columns = liftover_set.columns.str.replace("Macaca_mulatta", "macaque", regex=False)
liftover_set.columns = liftover_set.columns.str.replace("Callithrix_jacchus", "marmoset", regex=False)

## ortholog_set
ortholog_set = liftover_set[liftover_set["ortholog"] == True]
len(ortholog_set)

## Annos
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## Save for CACTUS analysis
liftover_set.to_csv(os.path.join(anno_dir, "human_ref_liftover_HALPER_minMatch_0-5.tsv"), sep="\t", index=False, header=True)