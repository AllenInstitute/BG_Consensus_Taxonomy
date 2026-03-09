import pandas as pd
import pyranges as pr
import anndata as ad
import os

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
from matplotlib.colors import ListedColormap

##
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Load variant library
HAR_lib = pd.read_excel(os.path.join(variant_dir, "HARs_312_hg38_Science2023.xlsx"), sheet_name="zooHARs")
HAR_lib.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "simple_name": "HAR_ID"}, inplace=True)
HAR_pr = pr.PyRanges(HAR_lib.loc[:, ["Chromosome", "Start", "End", "HAR_ID"]])

hCONDELs_lib = pd.read_excel(os.path.join(variant_dir, "hCONDELs_10032_2023Science.xlsx"), sheet_name="hCONDEL metadata")
hCONDELs_lib.rename(columns={"hg38_cons_chr": "Chromosome", "hg38_cons_start_pos": "Start", "hg38_cons_end_pos": "End"}, inplace=True)
hCONDELs_pr = pr.PyRanges(hCONDELs_lib.loc[:, ["Chromosome", "Start", "End", "hCONDEL_ID"]])

HAQER_lib = pd.read_csv(os.path.join(variant_dir, "haqer.hg38.bed"), sep="\t", header=None)
HAQER_lib.columns = ["Chromosome", "Start", "End", "HAQER_ID", "Length"]
HAQER_pr = pr.PyRanges(HAQER_lib)

## Load in peak x group boolean matrix
species = "human"
peak_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/{species}/ATAC"
peak_group = pd.read_csv(os.path.join(peak_dir, "Group_by_peaks.csv"), index_col=0)

##
peak_group["Chromosome"] = peak_group.index.str.split(":").str[0]
peak_group["Start"] = peak_group.index.str.split(":").str[1].str.split("-").str[0].astype(int)
peak_group["End"] = peak_group.index.str.split(":").str[1].str.split("-").str[1].astype(int)

## Fix names
peak_group.columns = [c.replace(" ", "_") for c in peak_group.columns]
peak_pr = pr.PyRanges(peak_group.loc[:, ["Chromosome", "Start", "End"] + list(peak_group.columns[:-3])])

## ----------------------------------
## HAR analysis
## ----------------------------------
overlap_HAR = peak_pr.join(HAR_pr, suffix="_HAR")
overlap_HAR_df = overlap_HAR.df
overlap_HAR_df.set_index("HAR_ID", inplace=True)

## ----------------------------------
## hCONDEL analysis
## ----------------------------------
overlap_hCONDELs = peak_pr.join(hCONDELs_pr, suffix="_hCONDEL")
overlap_hCONDELs_df = overlap_hCONDELs.df
overlap_hCONDELs_df.set_index("hCONDEL_ID", inplace=True)

## ----------------------------------
## HAQER analysis
## ----------------------------------
overlap_HAQER = peak_pr.join(HAQER_pr, suffix="_HAQER")
overlap_HAQER_df = overlap_HAQER.df
overlap_HAQER_df.set_index("HAQER_ID", inplace=True)

## Check for intersection with ABC scores
ABC_scores_f = pd.read_csv(os.path.join(analysis_dir, "atac", "abc-model", f"ABC_scores_{species}_filtered_ABC_gt_0.02.csv"))
ABC_scores_pr = pr.PyRanges(ABC_scores_f.rename(columns={"chrom": "Chromosome", "peak_start": "Start", "peak_end": "End"}).loc[:, ["Chromosome", "Start", "End", "gene_id", "ABC"]])

## 
overlap_HAQER_ABC = overlap_HAQER.join(ABC_scores_pr, suffix="_ABC").df
haqer_genes = overlap_HAQER_ABC.gene_id.unique()

## Human specialized genes (hDEGs)
de_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/RNA/"
human_specific_df = pd.read_csv(os.path.join(de_dir, "human_specific_DE_genes_bg.csv"))

##
human_specific_genes = set(human_specific_df.gene.unique())
overlap_HAQER_ABC["hDEG"] = overlap_HAQER_ABC.gene_id.isin(human_specific_genes)

## Summarize for plot
overlap_HAQER_df["overlaps_hDEG"] = overlap_HAQER_df.index.isin(overlap_HAQER_ABC.HAQER_ID[overlap_HAQER_ABC.hDEG])

import gseapy as gp

# prefetch gene sets to make faster
gene_set_library = 'GO_Biological_Process_2025'
gs_dict = gp.get_library(name=gene_set_library)  # dict: term -> [genes]
if not isinstance(gs_dict, dict): #in case it's a lol
    gs_dict = {k: v for k, v in gs_dict}
 
gs_dict = {k:v for k,v in gs_dict.items() if len(gs_dict[k])>5} #by hand filter
gs_dict = {k:v for k,v in gs_dict.items() if len(gs_dict[k])<100}

##
gene_list = overlap_HAQER_ABC.gene_id[overlap_HAQER_ABC.hDEG].unique().tolist()
enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=gs_dict,
                 organism='human',
                 # background=human_specific_df.gene.unique(),
                 outdir=None,
                )

res_df = enr.results.sort_values(["Adjusted P-value", "P-value"]).reset_index(drop=True)
res_sig_df = res_df.loc[res_df["Adjusted P-value"] < 0.05, :]


## ----------------------------------
## Load in annotations and data
## ----------------------------------
merged_adata = ad.read_h5ad(os.path.join(analysis_dir, "cactus", "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"))
merged_adata = merged_adata[merged_adata.obs["species"] == "human", :]

human_adata = ad.read_h5ad(os.path.join(data_dir, "human", "crested_adata", "human_basalganglia_hmba_pre-print_crested_500bp.h5ad"))

## ----------------------------------
## Plot individual heatmaps
## ----------------------------------

## params
add_expr = False

## Plot heatmaps
for df, variant_type, meta_cols in [
    (overlap_HAR_df, "HAR", ["Chromosome", "Start", "End", "Start_HAR", "End_HAR", "HAR_ID"]),
    (overlap_hCONDELs_df, "hCONDEL", ["Chromosome", "Start", "End", "Start_hCONDEL", "End_hCONDEL", "hCONDEL_ID"]),
    (overlap_HAQER_df, "HAQER", ["Chromosome", "Start", "End", "Start_HAQER", "End_HAQER", "HAQER_ID", "Length"]),
]:
    ## Filter to cell-type columns
    celltype_cols = [c for c in df.columns if c not in meta_cols]
    celltype_cols = [c for c in celltype_cols if c in human_adata.obs_names]
    ## Update df to help with bringing in ATAC data
    df["variant-id"] = df.index
    df["region-id"] = df["Chromosome"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    df.set_index("region-id", inplace=True)
    ## Create matrix
    mat = df[celltype_cols].astype(int)
    dup_idx = mat.index.duplicated()
    mat = mat.loc[~dup_idx, :]
    ## Sort rows by total activity
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index, 
                    cluster_meta.Group[cluster_meta.Group.isin(celltype_cols)].values] 
    activity_sum = mat.sum(axis=1)
    mat = mat.loc[activity_sum > 0, :]
    mat = mat.sort_values(mat.columns.tolist(), ascending=False)
    if add_expr:
        ## Add in accessibility values
        atac_values = human_adata[celltype_cols, mat.index].X.toarray().T  # peaks x cell types
        ## set all mat elements to 1
        mat = mat * atac_values
        ## Min-max normalize each peak
        mat = (mat - mat.min(axis=1).values.reshape(-1,1)) / ((mat.max(axis=1) - mat.min(axis=1)).values.reshape(-1,1))
        mat.fillna(0, inplace=True)
    ## Sort rows by number of active cell types
    mat["activity_sum"] = activity_sum[mat.index]
    ## Plot
    plt.figure(figsize=(8, 10))
    sns.heatmap(
        mat[cluster_meta.Group[cluster_meta.Group.isin(celltype_cols)].values].T,
        cmap="Greys" if not add_expr else "viridis",
        cbar=False,
        linewidths=0,
        linecolor=None,
        yticklabels=True,
        xticklabels=False,
    )
    plt.title(f"{variant_type} Overlaps across Cell Types", fontsize=14)
    plt.xlabel("Peaks")
    plt.ylabel("Cell Types")
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "atac", "har_haqer_hcondel", f"{variant_type}_expr_{add_expr}_overlap_heatmap.png"), dpi=900)
    ## Plot the Activity Sum as a heatmap
    plt.figure(figsize=(2, 10))
    sns.heatmap(
        mat[["activity_sum"]],
        cmap="Reds",
        cbar=True,
        linewidths=0,
        linecolor=None,
        yticklabels=False,
        xticklabels=False,
    )
    plt.title(f"{variant_type} Overlap Activity Sum", fontsize=14)
    plt.xlabel("Peaks")
    plt.ylabel("Activity Sum")
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "atac", "har_haqer_hcondel", f"{variant_type}_overlap_activity_sum_heatmap.png"), dpi=900)

## -------------------------------------------
## Plot specificty and accessibility values
## -------------------------------------------

variant_type = "HAQER"

##
haqer_gini = merged_adata.obs.loc[mat.index, "gini_scores"]
plt.figure(figsize=(2, 10))
sns.heatmap(
    haqer_gini.values.reshape(1, -1).T,
    cmap="Blues",
    cbar=True,
    linewidths=0,
    linecolor=None,
    yticklabels=False,
    xticklabels=False,
)
plt.title(f"{variant_type} Overlap Gini Score", fontsize=14)
plt.xlabel("Peaks")
plt.ylabel("Gini Score")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "atac", "har_haqer_hcondel", f"{variant_type}_overlap_gini_score_heatmap.png"), dpi=900)

##
haqer_N1 = merged_adata.obs.loc[mat.index, "evo_distance"]
plt.figure(figsize=(2, 10))
sns.heatmap(
    haqer_N1.values.reshape(1, -1).T,
    cmap="Purples",
    cbar=True,
    linewidths=0,
    linecolor=None,
    yticklabels=False,
    xticklabels=False,
)
plt.title(f"{variant_type} Overlap N1 Score", fontsize=14)
plt.xlabel("Peaks")
plt.ylabel("N1 Score")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "atac", "har_haqer_hcondel", f"{variant_type}_overlap_N1_score_heatmap.png"), dpi=900)

##
overlap_HAQER_df = overlap_HAQER_df.loc[mat.index, :]
cmap = ListedColormap(["lightgrey", "darkgreen"])
plt.figure(figsize=(2, 10))
sns.heatmap(
    overlap_HAQER_df["overlaps_hDEG"].values.reshape(1, -1).T,
    cmap=cmap,
    cbar=True,
    linewidths=0,
    linecolor=None,
    yticklabels=False,
    xticklabels=False,
)
plt.title(f"{variant_type} Overlap hDEG Score", fontsize=14)
plt.xlabel("Peaks")
plt.ylabel("hDEG Score")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "atac", "har_haqer_hcondel", f"{variant_type}_overlap_hDEG_score_heatmap.png"), dpi=900)
