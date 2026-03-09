import os
import pandas as pd
import numpy as np

## Handle nans
def parse_coordinates(coord):
    if pd.isna(coord) or not isinstance(coord, str) or ":" not in coord or "-" not in coord:
        return pd.Series([np.nan, np.nan, np.nan])
    try:
        chrom, pos = coord.split(":")
        start, end = pos.split("-")
        return pd.Series([chrom, int(start), int(end)])
    except Exception:
        return pd.Series([np.nan, np.nan, np.nan])

## -------------------------
##
## -------------------------
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"

## -------------------------
enhancer_table = pd.read_excel(os.path.join(analysis_dir, "crested", "enhancer_tables", "striatal_tools_supp_table1.xlsx"))
enhancer_table["enhancer ID"] = enhancer_table["enhancer ID"].str.strip()

## Coordinates
enhancer_table[["Chromosome", "Start", "End"]] = enhancer_table["Genome coordinates (cloned)"].apply(parse_coordinates)

## Remove NA
enhancer_table = enhancer_table.dropna(subset=["Chromosome", "Start", "End"]).reset_index(drop=True)

##
enhancer_bed = enhancer_table.loc[:, ["Chromosome", "Start", "End", "Fragment cloned size (bp)", "Species", "enhancer ID"]]
enhancer_bed.rename(columns={ "Fragment cloned size (bp)": "size", "Species": "species", "enhancer ID": "enhancer_id" }, inplace=True)

##
bed_output_path = os.path.join(analysis_dir, "crested", "enhancer_tables", "striatal_tools_enhancers.bed")
enhancer_bed.to_csv(bed_output_path, sep="\t", header=False, index=False)
