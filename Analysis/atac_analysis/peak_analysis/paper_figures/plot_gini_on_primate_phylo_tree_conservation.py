import os
from collections import defaultdict

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

from Bio import Phylo

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"
atac_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables"

## Load species and ancestral .h5ad
adata_files = {
    "species": ad.read_h5ad(os.path.join(analysis_dir, "CACTUS", "all_species_zoonomia_overlap_HALPER_NCBI_peak_categories_with_gini.h5ad")),
    "ancestral": ad.read_h5ad(os.path.join(analysis_dir, "CACTUS", "all_species_zoonomia_overlap_HALPER_ancestral.h5ad"))
}

## Merge all adatas
merged_adata = ad.concat(
    adata_files.values(),
    axis=1,           
    join="outer", ## include all features (vars)
    merge="unique",   
)
del merged_adata.obs['source']
merged_adata.obs_names_make_unique()

## Compute alignment 90%
peak_length = 501  # Fixed length from HALPER
merged_adata.layers[">=90%_aligned"] = (merged_adata.X/peak_length >= 0.9).astype(int)

## --------------------------------------------
## Setup tree for evo distance walking
## --------------------------------------------
tree = Phylo.read(os.path.join(cactus_path, "447-mammalian-2022v1.nh"), "newick")  # or .newick
tree.ladderize()  # optional tidy ordering

def tree_to_dataframe(tree):
    records = []
    ##
    def walk(clade, parent_name=None):
        node_name = clade.name if clade.name else f"internal_{len(records)}"
        is_leaf = len(clade.clades) == 0
        records.append({
            "node_name": node_name,
            "parent_name": parent_name,
            "branch_length": clade.branch_length,
            "is_leaf": is_leaf
        })
        for child in clade.clades:
            walk(child, node_name)
    ##
    walk(tree.root)
    return pd.DataFrame(records)

phylo_df = tree_to_dataframe(tree)
phylo_df["branch_length"] = phylo_df["branch_length"].fillna(0)

## Subset internal nodes only (if that's what alignment uses)
internal_df = phylo_df[~phylo_df["is_leaf"]].reset_index(drop=True)
internal_df["branch_length_norm"] = internal_df["branch_length"] / phylo_df["branch_length"].sum()

print(internal_df.head())
print(f"Internal nodes: {len(internal_df)} / Total nodes: {len(phylo_df)}")

## ---------------------------------------------------------------------
## Primate specific conservaion gini plotting
## ---------------------------------------------------------------------

primate_categories = pd.read_csv(os.path.join(anno_dir, "primates_CACTUS_447_with_ncbi_taxa_custom_names.csv"), index_col=0)

## Define species set for pre-defined collections
species_set = primate_categories.index.tolist()
mrca = tree.common_ancestor(species_set)

## Get all children, leaves and internal nodes of the primate MRCA
primate_leaves = [leaf.name for leaf in mrca.get_terminals()]
print(f"Primate leaves: {primate_leaves}")

primate_internal = [internal.name for internal in mrca.get_nonterminals()]
print(f"Primate internal nodes: {primate_internal}")

node_collection = primate_leaves + primate_internal
print(f"Total primate nodes: {len(node_collection)}")

## ---------------------------------------------------------------------
## Primate specific evo distance computation
## ---------------------------------------------------------------------

human_adata = merged_adata[merged_adata.obs["species"] == "human"] 

## Setup alignment and evo distances (branch lengths)
alignment_matrix = human_adata[:,primate_internal].layers["alignment_90pct"].toarray()  # peaks × nodes

primate_df = internal_df[internal_df["node_name"].isin(primate_internal)].set_index("node_name")
branch_lengths = primate_df.loc[primate_internal, "branch_length_norm"].to_numpy()

## Compute evo distances
human_adata.obs["evo_distance"] = (alignment_matrix * branch_lengths).sum(axis=1)

## ---------------------------------------------------------------------
## Plot gini vs evo distance for primate conserved peaks
## ---------------------------------------------------------------------

gini_scores = human_adata.obs["gini_scores"]
evo_distances = human_adata.obs["evo_distance"]

# plt.figure(figsize=(6,4))
# plt.scatter(evo_distances, gini_scores, s=0.001, alpha=0.5, color="#4dac26")
# plt.xlabel("Primate Conservation Evo Distance")
# plt.ylabel("Minimum Gini Score")
# plt.title("Primate Conservation vs Gini Score")
# plt.xlim(0, evo_distances.max()*1.05)
# plt.ylim(0, 1.0)
# plt.axhline(y=0.7, color='black', linestyle='--')
# plt.axvline(x=0.25, color='black', linestyle='--')
# plt.tight_layout()
# plt.savefig(os.path.join(analysis_dir, "specificity", "figures", "primate_conservation_gini_scatter.png"), dpi=600)
# plt.close()

manual_set =    [
                    'Great_apes', 
                    'Small_apes', 
                    'Old_World_monkeys', 
                    'New_World_monkeys',
                    'Lemuriforms',
                    'Tarsiers'
                ] + \
                [
                    "Great_apes;Small_apes", 
                    "Old_World_monkeys;Great_apes;Small_apes", 
                    "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys",
                    "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers",
                    "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers;Lemuriforms"
                ] + \
                ["None"]

plot_adata = human_adata[human_adata.obs["ancestral_category"].isin(manual_set), :].copy()

summary_df = (
    plot_adata.obs
    .groupby("ancestral_category")["gini_scores"]
    .describe()  # gives count, mean, std, min, quartiles, max
    .sort_values("mean", ascending=False)
)
print(summary_df)


import seaborn as sns
plt.figure(figsize=(7, 4))
sns.boxplot(
    data=plot_adata.obs,
    x="ancestral_category",
    y="gini_scores",
    order=manual_set,  # optional ordering
)
plt.xticks(rotation=45, ha='right')
plt.ylabel("Gini score")
plt.xlabel("Primate clade")
plt.title("Distribution of Gini scores by primate clade")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "specificity", "figures", "primate_conservation_gini_boxplot.png"), dpi=600)



## Load in marker peak information
marker_peaks = pd.read_csv(os.path.join(analysis_dir, "specificity", "top_100_peaks_per_cell_type_per_species.csv"))

## Add cell type info to human_adata
human_adata.obs["cell_type"] = "Unknown"
for cell_type in marker_peaks["cell_type"].unique():
    peaks_for_type = marker_peaks[marker_peaks["cell_type"] == cell_type]["region"].tolist()
    human_adata.obs.loc[human_adata.obs_names.isin(peaks_for_type), "cell_type"] = cell_type

plot_adata = human_adata[human_adata.obs["ancestral_category"].isin(manual_set), :].copy()
plot_adata = plot_adata[(plot_adata.obs["cell_type"] != "Unknown"), :].copy()
plot_adata = plot_adata[plot_adata.obs["cell_type"].str.contains("Gluta-Dorsal"), :].copy()

plot_adata.obs["ancestral_category"] = pd.Categorical(
    plot_adata.obs["ancestral_category"],
    categories=manual_set,
    ordered=True
)

# Sort by cell_type then by assigned_clade order
plot_df = (
    plot_adata.obs
    .sort_values(["cell_type", "ancestral_category"])
)

plot_df = plot_df.loc[plot_df["gini_scores"] >0.75, :]


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
sns.violinplot(
    data=plot_df,
    x="ancestral_category",
    y="gini_scores",
    hue="cell_type",
    order=manual_set,
    split=True,          # 🟢 shows both cell types side-by-side in one violin
    inner="quartile",
    cut=0
)
plt.xticks(rotation=45, ha="right", fontsize=2)
plt.xlabel("Primate clade")
plt.ylabel("Gini score")
plt.title("Cell-type–specific Gini distributions across primate clades")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "specificity", "figures", "primate_conservation_gini_lines.png"), dpi=600)
