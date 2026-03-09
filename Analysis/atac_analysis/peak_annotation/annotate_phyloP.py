import os
import pandas as pd
import pyranges as pr
import pyBigWig
from tqdm import tqdm
import numpy as np
from tqdm import tqdm

## Import progress_apply
tqdm.pandas()

## Peak data
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
species_peaks = {
    "Homo_sapiens": pd.read_csv(os.path.join(data_dir, "human/ATAC/merged_peaks.bed"), 
                                    names=["Chromosome", "Start", "End", "ID"], sep="\t", header=None)
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

## -------------------------------
## Directories
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## PhyloP bigWigs (example paths, update to your downloads)
human_phyloP_file = os.path.join(work_dir, "anno_tables", "ATAC", "phyloP", "hg38.phyloP100way.bw")

## -------------------------------
## Helper: get mean phyloP per interval
def get_mean_phyloP(bw, chrom, start, end):
    try:
        vals = bw.values(chrom, int(start), int(end), numpy=True)
        if vals is None:
            return np.nan
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0:
            return np.nan
        return float(vals.mean())
    except RuntimeError:
        # Chromosome naming mismatch etc.
        return np.nan

## -------------------------------
## Main loop
species_peak_files = {**species_peaks, **species_liftovers}

## Open phyloP bigWig
with pyBigWig.open(human_phyloP_file) as bw:
    ##
    for species, peaks in species_peak_files.items():
        print(f"Processing {species}...")
        ## Annotate mean phyloP per peak
        tqdm.pandas(desc=f"{species} peaks")
        peaks["phyloP_mean"] = peaks.progress_apply(
            lambda r: get_mean_phyloP(bw, r["Chromosome"], r["Start"], r["End"]),
            axis=1
        )
        ## Save
        out_file = os.path.join(anno_dir, f"{species}_peaks_annotated_phyloP.csv")
        peaks.to_csv(out_file, index=False)
