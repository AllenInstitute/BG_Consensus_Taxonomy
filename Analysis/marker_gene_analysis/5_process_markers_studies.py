import pandas as pd
import numpy as np
import anndata as ad
import os, sys, glob, re

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript"
marker_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Accuracy_analysis"

## --------------------------------
## Now read in our markers sets
marker_files = glob.glob(os.path.join(marker_dir, "*_markers_092325.tsv"))

marker_sets = {}
for file in marker_files:
    species = os.path.basename(file).replace("_markers_092325.tsv", "")
    markers = []
    with open(file, "r") as f:
        for line in f:
            markers.append(line.strip())
    marker_sets[species] = markers
    print(f"{species}: {len(markers)} markers")

## --------------------------------
## Read in Yadav and Russ markers
yadav_markers = pd.read_excel(os.path.join(marker_dir, "markers_yadav.xlsx"))
russ_markers = pd.read_excel(os.path.join(marker_dir, "markers_russ.xlsx"))

## --------------------------------
## Proces wilcox derived markers
russ_markers_wilcox = russ_markers.iloc[0:20,:]
russ_markers_wilcox = russ_markers_wilcox.drop(columns=["WILCOX"])

## Collapse all markers across columns into a single list
russ_markers_list = []
for ct in russ_markers_wilcox.columns:
    print(f"{ct}: {russ_markers_wilcox[ct].nunique()} unique markers")
    markers = russ_markers_wilcox[ct].dropna().unique().tolist()
    russ_markers_list.extend(markers)

## Ensure markers are unique
marker_sets["russ_et_al"] = list(set(russ_markers_list))

## --------------------------------
## Process Yadav markers
yadav_markers_list = []
for ct in yadav_markers.columns:
    print(f"{ct}: {yadav_markers[ct].nunique()} unique markers")
    markers = yadav_markers[ct].dropna().unique().tolist()
    yadav_markers_list.extend(markers)

## Ensure markers are unique
marker_sets["yadav_et_al"] = list(set(yadav_markers_list))

## --------------------------------
## Save out processed markers

markers = pd.DataFrame({k: pd.Series(v) for k, v in marker_sets.items()})
markers = markers.dropna(how="all")

## Ensure all columns are capitalized
markers.columns = [col.capitalize() for col in markers.columns]
## Write out
markers.to_csv(os.path.join(marker_dir, "processed_marker_sets_092325.csv"))