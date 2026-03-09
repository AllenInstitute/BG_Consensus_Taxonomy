import pandas as pd

## Load in gini scores
gini_scores = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/specificity/gini_scores_combined.csv")

## Load in accessolog scores
accessolog_scores = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Accessologs/accessologs_mammalian_orthologs.csv")
accessolog_scores = accessolog_scores.drop_duplicates(subset=["region"])

## Merge data
merged = pd.merge(accessolog_scores, gini_scores, left_on="region", right_on="region", how="inner")

##