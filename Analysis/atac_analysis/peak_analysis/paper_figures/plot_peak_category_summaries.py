import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/ATAC-conservation/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"

## spinalcord metadata
cluster_meta = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

##
merged_adata = ad.read_h5ad(os.path.join(analysis_dir, "CACTUS", "all_species_zoonomia_overlap_HALPER_NCBI_peak_categories.h5ad"))

## Sequence categories
seq_cats = ['species-specific', 'seq-conserved', 'epi-conserved', 'epi-conserved-markers', 'human-biased']

## --------------------------------
## Plot promoter proportion for each sequence category
## --------------------------------
for seq_cat in seq_cats:
    plot_adata = merged_adata[(merged_adata.obs[seq_cat] == True), :].copy()
    plot_adata.obs["promoter"] = plot_adata.obs["promoter"].astype(str)
    ## Calculate promoter proportion
    promoter_prop = plot_adata.obs.groupby("promoter").size()
    promoter_prop = promoter_prop / promoter_prop.sum()
    ## Pie chart
    plt.figure(figsize=(4,4))
    plt.pie(promoter_prop, 
            colors=["olive", "blue"])
    plt.title(f"{seq_cat} (n={plot_adata.n_obs})")
    plt.savefig(os.path.join(analysis_dir, "peak_categories", f"{seq_cat}_promoter_proportion_piechart.png"))
    plt.close()

## --------------------------------
## Density plot colored by species of gini_scores
## --------------------------------
for seq_cat in seq_cats:
    plot_adata = merged_adata[(merged_adata.obs[seq_cat] == True), :].copy()
    plt.figure(figsize=(6,4))
    sns.kdeplot(data=plot_adata.obs, 
                x="gini_scores", 
                color="olive",
                common_norm=False, 
                fill=True, 
                alpha=0.5)
    plt.title(f"{seq_cat} (n={plot_adata.n_obs})")
    plt.xlabel("Gini Scores")
    plt.ylabel("Density")
    plt.savefig(os.path.join(analysis_dir, "peak_categories", f"{seq_cat}_gini_scores_density.png"))
    plt.close()

## --------------------------------
## Stacked barplot of TE elements for each sequence category
## --------------------------------
te_order = ['LTR', "LINE", "SINE", "DNA"]
for seq_cat in seq_cats:
    plot_adata = merged_adata[(merged_adata.obs[seq_cat] == True) & 
                                (merged_adata.obs["species"] == "human"), :].copy()
    te_counts = plot_adata.obs['repClass'].value_counts().reindex(te_order).fillna(0)
    te_props = te_counts / te_counts.sum()
    ## Barplot
    plt.figure(figsize=(6,4))
    te_props.plot(kind='bar', color=cm.tab20.colors)
    plt.title(f"{seq_cat} (n={plot_adata.n_obs})")
    plt.ylabel("Proportion")
    plt.xlabel("TE Class")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "peak_categories", f"{seq_cat}_TE_class_barplot.png"))
    plt.close()


te_order = ['LTR', "LINE", "SINE", "DNA"]
te_colors = [cm.tab20.colors[i] for i in range(len(te_order))] + [(1.0, 1.0, 1.0)]  # gray for "Other"

# Collect proportions for each seq_cat
te_prop_df = pd.DataFrame(index=te_order + ["Other"])

for seq_cat in ['species-specific', 'human-biased', 'seq-conserved', 'epi-conserved']:
    plot_adata = merged_adata[
        (merged_adata.obs[seq_cat] == True) &
        (merged_adata.obs["species"] == "human"),
        :
    ].copy()
    # Get counts for all TE classes
    te_counts_all = plot_adata.obs['repClass'].value_counts()
    # Split into main and other
    te_counts = te_counts_all.reindex(te_order).fillna(0)
    other_sum = te_counts_all[~te_counts_all.index.isin(te_order)].sum()
    te_counts["Other"] = other_sum
    # Convert to proportions
    te_props = te_counts / te_counts.sum()
    te_prop_df[seq_cat] = te_props

# Transpose for plotting
te_prop_df = te_prop_df.T

# Plot a single stacked bar plot
plt.figure(figsize=(6, 4))
bottom = np.zeros(len(te_prop_df))
for i, te in enumerate(te_prop_df.columns):
    plt.bar(te_prop_df.index, 
            te_prop_df[te], 
            bottom=bottom, 
            color=te_colors[i], 
            label=te)
    bottom += te_prop_df[te].values

plt.ylabel("Proportion")
plt.xlabel("Sequence category")
plt.xticks(rotation=45, ha='right')
plt.ylim(0, 0.6)  # <-- limit y-axis to 0.8 for more dynamic range
plt.legend(title="TE Class", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig(os.path.join(analysis_dir, "peak_categories", "all_seqcat_TE_class_stacked_bar.png"))
plt.close()