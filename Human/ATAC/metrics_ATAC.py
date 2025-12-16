import os, sys, re, glob
from tqdm import tqdm
tqdm.pandas()
## Model
import crested
import keras
import tensorflow as tf
## General data
import pandas as pd
import anndata as ad
import numpy as np
## Genomics handling
import pyranges as pr
import pysam

## Utility functions
sys.path.append("/code/")
from utils import *

## -------------------------
## Arguments
## -------------------------

anno_level = "class"

bigwigs_path = os.path.join("/data", "wb-aav-toolbox", "AIT21", anno_level, "bigwig")
peak_path  = os.path.join("/data", "wb-aav-toolbox", "AIT21", anno_level, "peaks", f"{anno_level}_peaks.bed")

## -------------------------
## Load prioritization table
## -------------------------

enhancer_table_path = f"/data/wb-aav-toolbox-priortized-peaks/AIT21/peakrankr/{anno_level}/{anno_level}_differential_peaks_ranked.bed"
enhancer_table = pd.read_csv(enhancer_table_path, sep="\t")

## -------------------------
## Load data
## -------------------------
adata = crested.import_bigwigs(
    bigwigs_folder=bigwigs_path,
    regions_file=peak_path, ## Ideally this would be peaks from the prioritization table
    target_region_width=None, ## Use peak sizes if variable.
    target="mean",
)
adata.X = np.nan_to_num(adata.X, nan=0.0)

## ------------------------
## Compute ATAC metrics
## ------------------------

adata = calc_gini(adata)

## Find index of max value per row and gather cell type
max_idx = np.argmax(adata.layers["gini_scores"], axis=0)
max_cell_types = adata.obs_names[max_idx]

## Get the max values
max_values = np.max(adata.layers["gini_scores"], axis=0)

## Build metric dataframe
metric_gini_specificity = pd.DataFrame({
    "max_celltype": max_cell_types,
    "max_value": max_values
}, index=adata.var_names)
metric_gini_specificity.head()

## Check overlap of enhancer_table and metric_gini_specificity index
overlap = enhancer_table.index.intersection(metric_gini_specificity.index)
print(f"Overlap between enhancer_table and metric_gini_specificity: {overlap}")

## Map results back to enhancer_table
enhancer_table = enhancer_table.join(metric_gini_specificity, how="left", rsuffix="_gini")