import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from Bio import Phylo

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/ATAC-conservation/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"

## spinalcord metadata
cluster_meta = pd.read_excel("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

##
merged_adata = ad.read_h5ad(os.path.join(analysis_dir, "CACTUS", "all_species_zoonomia_overlap_HALPER_NCBI_peak_categories_with_gini.h5ad"))

# ## primate ancestral
# ancestral_anno = pd.read_csv(os.path.join(analysis_dir, "CACTUS", "zoonomia_overlap_HALPER_ancestral_categories.csv"), index_col=0)

# ## Manual clade sets
# manual_set =    ["None"] + \
#                 [
#                     'Great_apes', 
#                     # 'Small_apes', 
#                     # 'Old_World_monkeys', 
#                     # 'New_World_monkeys',
#                     # 'Lemuriforms',
#                     # 'Tarsiers'
#                 ] + \
#                 [
#                     "Great_apes;Small_apes", 
#                     "Old_World_monkeys;Great_apes;Small_apes", 
#                     "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys",
#                     "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers",
#                     "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers;Lemuriforms"
#                 ]

# merged_adata.obs = merged_adata.obs.join(ancestral_anno, rsuffix="_cat", how="left")

# merged_adata.obs["ancestral_category"] = ancestral_anno.loc[merged_adata.obs_names, "ancestral_category"]
# merged_adata.obs["ancestral_category"][merged_adata.obs["ancestral_category"].isna()] = "None"

# ## Order sets by evo distance
# ordered_by_evo_distance = [
#     'Primates',        # humans
#     'Great_apes', 
#     'Small_apes', 
#     'Old_World_monkeys', 
#     'New_World_monkeys',
#     'Lemuriforms',
#     'Tarsiers'
#     'Dermoptera',      # colugos (closest living relatives to primates)
#     'Scandentia',      # treeshrews (next closest)
#     'Rodentia',        # rodents
#     'Lagomorpha',      # rabbits
#     'Eulipotyphla',    # shrews, hedgehogs
#     'Chiroptera',      # bats
#     'Carnivora',       # dogs, cats, seals
#     'Pholidota',       # pangolins (sister to Carnivora)
#     'Perissodactyla',  # horses, rhinos
#     'Artiodactyla',    # cows, pigs, whales
#     'Macroscelidea',   # elephant shrews
#     'Hyracoidea',      # hyraxes
#     'Proboscidea',     # elephants
#     'Sirenia',         # manatees, dugongs
#     'Cingulata',       # armadillos
#     'Pilosa',          # sloths, anteaters
#     'Tubulidentata'    # aardvark (most basal of these)
# ]

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

## Find human leaf
human_leaf = [term for term in tree.get_terminals() if "Homo_sapiens" in term.name]
human_leaf = human_leaf[0]

## Compute distance from human for each node
node_distances = {}
for clade in tree.find_clades(order='preorder'):
    node_name = clade.name if clade.name else f"internal_{id(clade)}"
    node_distances[node_name] = tree.distance(clade, human_leaf)

## Add to dataframe
phylo_df["distance_from_human"] = phylo_df["node_name"].map(node_distances)

## Sequence categories
seq_cats = ['species-specific', 'seq-conserved', 'epi-conserved', 'epi-conserved-markers', 'human-biased']

## --------------------------------
## Box and whisker of gini scores across manual sets
## --------------------------------
# gini_df = pd.DataFrame(columns=["gini", "species", "manual_set"])
# for species in merged_adata.obs["species"].unique():
#     for mset in ordered_by_evo_distance:
#         if mset not in merged_adata.obs.columns:
#             continue
#         plot_adata = merged_adata[(merged_adata.obs["species"] == species) & 
#                                     (merged_adata.obs[mset] == True), :].copy()
#         gini_scores = plot_adata.obs["gini_scores"]
#         gini_df = pd.concat([gini_df, pd.DataFrame({
#             "gini": gini_scores,
#             "species": [species] * len(gini_scores),
#             "manual_set": [mset] * len(gini_scores)
#         })], ignore_index=True)

for species in tqdm(merged_adata.var_names.unique()):
    merged_adata.obs[f"{species}_highly_conserved"] = ((merged_adata[:,species].X.todense() / 501) >= 1.0)

gini_df = pd.DataFrame(columns=["gini", "species", "manual_set"])
for ref_species in ["human"]:
    for species in tqdm(merged_adata.var_names.unique()):
        plot_adata = merged_adata[(merged_adata.obs["species"] == ref_species) & 
                            (merged_adata.obs[f"{species}_highly_conserved"] == True), :].copy()
        median_gini = plot_adata.obs.loc[:, "gini_scores"].median()
        print(f"{ref_species} - {species}: {median_gini}")
        gini_df = pd.concat([gini_df, pd.DataFrame({
            "gini": median_gini,
            "ref_species": ref_species,
            "query_species": species, 
        }, index=[0])], ignore_index=True)

##
# Merge gini_df with phylo_df to include distance info
plot_df = gini_df.merge(
    phylo_df[["node_name", "distance_from_human"]],
    left_on="query_species",
    right_on="node_name",
    how="left"
)

# Sort by distance
plot_df = plot_df.sort_values("distance_from_human")

# Plot
plt.figure(figsize=(10,4))
sns.lineplot(
    data=plot_df,
    x="query_species",
    y="gini",
    sort=False,
    marker="o"
)

plt.xticks(rotation=90)
plt.xlabel("Species (ordered by distance from human)")
plt.ylabel("Median Gini score")
plt.xticks([], [])
plt.title("Gini vs. Evolutionary Distance from Human")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "peak_categories", "gini_vs_evolutionary_distance_from_human.png"), dpi=300)





# # Filter to human or include all species as desired
# # plot_df = gini_df.loc[(gini_df.species == "human") & (gini_df.gini > 0.6), :].copy()
# # plot_df["manual_set"] = pd.Categorical(plot_df["manual_set"], categories=ordered_by_evo_distance, ordered=True)
# plot_df["query_species"] = pd.Categorical(plot_df["query_species"], categories=merged_adata.var_names.unique(), ordered=True)

# plt.figure(figsize=(10, 8))

# # Base boxplot
# sns.boxplot(
#     data=plot_df,
#     x="manual_set",
#     y="gini",
#     hue="species",
#     showcaps=True,
#     fliersize=0,       # hide default fliers
#     boxprops={'zorder': 2},
#     linewidth=1
# )

# # Add jittered points (show outliers and distribution)
# sns.stripplot(
#     data=plot_df,
#     x="manual_set",
#     y="gini",
#     hue="species",
#     dodge=True,        # separate by hue
#     alpha=0.4,
#     jitter=0.25,
#     zorder=1,
#     size=3
# )

# # Adjust legend to avoid duplication
# handles, labels = plt.gca().get_legend_handles_labels()
# plt.legend(handles[:len(plot_df['species'].unique())], labels[:len(plot_df['species'].unique())],
#            title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.xticks(rotation=45, ha='right')
# plt.ylabel("Gini score")
# plt.xlabel("Ancestral category")
# plt.title("Gini scores across ancestral categories (with jittered outliers)")
# plt.tight_layout()

# plt.savefig(os.path.join(analysis_dir, "peak_categories", "gini_scores_across_ancestral_categories_boxplot_jitter.png"))
# plt.close()