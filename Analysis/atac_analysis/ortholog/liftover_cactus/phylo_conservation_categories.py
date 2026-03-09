## --------------------------------------------------
## Plot phylogenetic tree of species in cactus HAL
## --------------------------------------------------
import os
import pandas as pd
import anndata as ad

import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

from Bio import Phylo

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
# cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Load species and ancestral .h5ad
adata_files = {
    "species": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad")),
    "ancestral": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER_ancestral.h5ad"))
}

## Merge all adatas
merged_adata = ad.concat(
    adata_files.values(),
    axis=1,           
    join="outer", ## include all features (vars)
    merge="unique",   
)
merged_adata.obs_names_make_unique()

## Compute alignment 90%
peak_length = 501  # Fixed length from HALPER
merged_adata.layers[">=90%_aligned"] = (merged_adata.X/peak_length >= 0.9).astype(int)

## --------------------------------------------
## Categorization of peaks based on pre-defined ancestral nodes
## --------------------------------------------
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables"
primate_categories = pd.read_csv(os.path.join(anno_dir, "primates_CACTUS_447_with_ncbi_taxa_custom_names.csv"), index_col=0)

tree = Phylo.read(os.path.join(cactus_path, "zoonomia_447_no_length.nh"), "newick")  # or .newick

## Define species set for pre-defined collections
primate_mrca = {}
for primates_shorthand in primate_categories["primates_shorthand"].unique():
    species_set = primate_categories.loc[primate_categories.primates_shorthand == primates_shorthand].index.tolist()
    ## Get the MRCA node
    mrca = tree.common_ancestor(species_set)
    print(f"{primates_shorthand}: {mrca.name}")
    ##
    primate_mrca[primates_shorthand] = mrca.name

order_mrca = {}
for order in adata_files["species"].var["order"].unique():
    species_set = adata_files["species"].var_names[adata_files["species"].var["order"] == order].tolist()
    ## Get the MRCA node
    mrca = tree.common_ancestor(species_set)
    print(f"{order}: {mrca.name}")
    ##
    order_mrca[order] = mrca.name

##
all_mrca = {**primate_mrca, **order_mrca}
all_mrca = {k: v for k, v in all_mrca.items() if pd.notna(k) and pd.notna(v)}

## Figure out which peaks align to the pre-defined mrca nodes
for category, mrca_node in all_mrca.items():
    if mrca_node in merged_adata.var_names:
        aligned = merged_adata[:, mrca_node].layers[">=90%_aligned"].toarray().flatten()
        merged_adata.obs[f"{category}"] = aligned
    else:
        print(f"MRCA node {mrca_node} for category {category} not found in adata.var_names")

## For each peak summarize its categories
def categorize_peak(row, sets=all_mrca):
    categories = [cat for cat in sets.keys() if row[cat]]
    if categories:
        return ";".join(categories)
    else:
        return "None"

## Apply categorization  
merged_adata.obs["ancestral_primate_category"] = merged_adata.obs.apply(categorize_peak, sets=primate_mrca, axis=1)
merged_adata.obs["ancestral_category"] = merged_adata.obs.apply(categorize_peak, axis=1)

# ## Write out
# merged_adata.obs.to_csv(os.path.join(analysis_dir, "analysis", "zoonomia_overlap_HALPER_ancestral_categories.csv"), index=True)

## --------------------------------------------
## Some plots
## --------------------------------------------

## Extract all known categories
categories = list(all_mrca.keys())

## Create boolean indicators for each category
for cat in categories:
    merged_adata.obs[cat] = merged_adata.obs["ancestral_category"].str.contains(cat, regex=False)

## Filter merged_adata to only categories with > 1000 peaks
category_counts = merged_adata.obs["ancestral_category"].value_counts()
filtered_categories = category_counts[category_counts > 1000].index.tolist()

##
for species_plot in ["human", "macaque", "marmoset"]:
    plot_adata = merged_adata[merged_adata.obs["ancestral_category"].isin(filtered_categories) & (merged_adata.obs["species"] == species_plot)].copy()
    ## Plot counts
    plt.figure(figsize=(10,6))
    sns.countplot(
        data=plot_adata.obs,
        y="ancestral_category",
        order=plot_adata.obs["ancestral_category"].value_counts().loc[filtered_categories].index,
        palette="viridis",
        stat="probability"
    )
    plt.xlabel("Number of Peaks")
    plt.ylabel("Ancestral Category")
    plt.title(f"Peak counts by ancestral category in {species_plot.capitalize()}")
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "figures", f"all_mrca_bar_{species_plot}.png"), dpi=900)
    plt.close()

##
for species_plot in ["human", "macaque", "marmoset"]:
    plot_adata = merged_adata[merged_adata.obs["ancestral_category"].isin(filtered_categories) & (merged_adata.obs["species"] == species_plot)].copy()
    ## Build the boolean DataFrame
    bool_df = plot_adata.obs[categories].astype(bool)
    upset_data = from_indicators(bool_df.columns, bool_df)
    ## Plot UpSet
    plt.figure(figsize=(8,6))
    UpSet(upset_data, subset_size='count', show_counts=True, sort_by='degree').plot()
    plt.suptitle("Peak overlap across primate MRCA categories")
    plt.savefig(os.path.join(analysis_dir, "figures", f"primate_mrca_upset_{species_plot}.png"), dpi=900)


## Manual clade sets
manual_set =  ["Great_apes",
                "Great_apes;Small_apes", 
                "Old_World_monkeys;Great_apes;Small_apes", 
                "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys",
                "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers",
                "Old_World_monkeys;Great_apes;Small_apes;New_World_monkeys;Tarsiers;Lemuriforms"] + ["None"]


##
for species_plot in ["human", "macaque", "marmoset"]:
    plot_adata = merged_adata[merged_adata.obs["ancestral_category"].isin(manual_set) & (merged_adata.obs["species"] == species_plot)].copy()
    print(plot_adata.obs["ancestral_category"].value_counts())
    ## Build the boolean DataFrame
    bool_df = plot_adata.obs[categories].astype(bool)
    upset_data = from_indicators(bool_df.columns, bool_df)
    ## Plot UpSet
    plt.figure(figsize=(10,6))
    UpSet(upset_data, subset_size='count', show_percentages="{:.1%}", sort_by='degree').plot()
    plt.suptitle("Peak overlap across primate MRCA categories")
    plt.savefig(os.path.join(analysis_dir, "figures", f"primate_mrca_upset_{species_plot}_manual_cat.pdf"), dpi=900)

## human
species_plot = "human"
plot_adata = merged_adata[merged_adata.obs["ancestral_category"].isin(manual_set) & (merged_adata.obs["species"] == species_plot)].copy()
print(plot_adata.obs["ancestral_category"].value_counts())

## Build the boolean DataFrame
bool_df = pd.read_csv("bool_df.csv", index_col=0) #plot_adata.obs[categories].astype(bool)
upset_data = from_indicators(bool_df.columns, bool_df)

## Load in peak x group boolean matrix
species = "human"
peak_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/{species}/ATAC"
peak_group = pd.read_csv(os.path.join(peak_dir, "Group_by_peaks.csv"), index_col=0)

## Subset peak group by plot_adata
peak_subset = peak_group.loc[bool_df.index]

## Assign each peak to either Neuron, Non-neuron, or both
def assign_celltype(row):
    neuron_groups = cluster_meta.loc[cluster_meta.Neighborhood != "Nonneuron", "Group"].tolist()
    neuron_groups = [g for g in neuron_groups if g in peak_subset.columns]
    non_neuron_groups = cluster_meta.loc[cluster_meta.Neighborhood == "Nonneuron", "Group"].tolist()
    non_neuron_groups = [g for g in non_neuron_groups if g in peak_subset.columns]
    is_neuron = row[neuron_groups].any()
    is_non_neuron = row[non_neuron_groups].any()
    ##
    if is_neuron and is_non_neuron:
        return "Both"
    elif is_neuron:
        return "Neuron"
    elif is_non_neuron:
        return "Non-neuron"
    else:
        return "Unassigned"
    
##
peak_subset["celltype"] = peak_subset.apply(assign_celltype, axis=1)

# Add celltype annotation (categorical)
upset_data["celltype"] = peak_subset["celltype"].values

# Define a consistent palette
palette = sns.color_palette("tab20", n_colors=peak_subset["celltype"].nunique())
palette_dict = dict(zip(peak_subset["celltype"], palette))

# ---- Aggregate by UpSet index & cell type ----
# Count how many of each cell type occur in each subset
composition_df = (
    upset_data.groupby(level=list(range(upset_data.index.nlevels)))["celltype"]
    .value_counts()
    .rename("count")
    .reset_index()
)

## Plot UpSet
upset = UpSet(
    upset_data,
    sort_by="degree",
    show_counts=False,
    show_percentages=None,
    intersection_plot_elements=0,
)
colors = ["#91f4bb", "#19613b", "#a8afa5"]
# Replace top bars with stacked composition bars
upset.add_stacked_bars(
    by="celltype", colors=colors,  title="Count by cell type"
)

# ---- Final plot ----
plt.figure(figsize=(10, 6))
upset.plot()
plt.suptitle(f"Cell-type composition of peak overlap subsets")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "figures", f"primate_mrca_upset_celltype_comp.pdf"), dpi=900)