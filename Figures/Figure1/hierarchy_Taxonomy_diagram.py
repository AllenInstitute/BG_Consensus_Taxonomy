import sys, os
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py

## Helpful locations which are assumed to already exist
annodir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
datadir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
# manuscript_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/figure1"

## -------------------------
## AnnoTable
## -------------------------

##### -- Load Data -- ######
# cross_species_paths = {
#     "HMBA:Human" : os.path.join(datadir, "human", "Human_HMBA_basalganglia_AIT_pre-print.h5ad"),
#     "HMBA:Macaque" : os.path.join(datadir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print.h5ad"),
#     "HMBA:Marmoset": os.path.join(datadir, "marmoset", "Marmoset_HMBA_basalganglia_AIT_pre-print.h5ad"),
# }

# species_obs = {}
# for species, path in cross_species_paths.items():
#     print(f"Processing {species}...")
#     with h5py.File(path) as f:
#         species_obs[species] = ad.experimental.read_elem(f['obs'])

# adata_obs = pd.concat(species_obs, axis=0)
# adata_obs = adata_obs.reset_index(level=0, drop=True)
# adata_obs.to_csv(os.path.join(datadir, "xspecies", "consensus_hmba_basalganglia_metadata.csv"))
adata_obs = pd.read_csv(os.path.join(datadir, "xspecies", "consensus_hmba_basalganglia_metadata.csv"), index_col=0)

anno = pd.read_excel(os.path.join(annodir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
anno = anno.drop_duplicates(subset=["Group"]).copy()
anno.reset_index(drop=True, inplace=True)
anno = anno.iloc[anno.display_order_group-1,:]

colors = pd.read_excel(os.path.join(annodir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="value_set_colors")

##### -- Order -- #####
## This order was determined first by a dendrogram using the top 1K marker genes then tinkered with manualy to group CHAT+ types.
ct_order = anno.Group.tolist()
anno_level = "Group"

################--------------- (1-1) region distribution pattern across Groups in Human ----------------################

anno_order_color = dict(zip(colors.label, colors.color_hex_triplet))

## Adjust roi names
adata_obs['anatomical_region_plot'] = adata_obs['anatomical_region']
adata_obs['anatomical_region_plot'] = adata_obs['anatomical_region_plot'].replace({
    'Br': 'Brain',
    'GPeC': 'GPe',
    'GPeR': 'GPe',
    'SN-VTA': 'SN',
})

anno_plot = "anatomical_region_plot"
# for species in adata_obs.organism.unique():
# print("Processing: " + species)
## Filter to species
species_obs = adata_obs[adata_obs['organism'] != 'Marmoset'] ## Remove marmoset as the tiles don't have ROI annotations
anno_table = species_obs[[anno_level, anno_plot]]
## Remove unused terms
anno_table = anno_table.loc[anno_table[anno_plot].isin(anno_order_color.keys()), :].copy()
## Define the desired order for 'anatomical_region'
anno_table[anno_plot] = pd.Categorical(anno_table[anno_plot], categories=anno_order_color.keys(), ordered=True)
## Define the desired order for anno_level
anno_table[anno_level] = anno_table[anno_level].astype(str).copy()
anno_table[anno_level] = pd.Categorical(anno_table[anno_level], categories=ct_order, ordered=True)
## Calculate the distribution of 'dissection' within each anno_level
distribution = pd.crosstab(anno_table[anno_level], anno_table[anno_plot])
## Reindex to include all groups, even those with no data
distribution = distribution.reindex(index=ct_order, fill_value=0)
## Convert counts to percentages
distribution_percentage = distribution.div(distribution.sum(axis=1), axis=0) * 100
## Gather only segments and colors for this species
colors_for_species = custom_colors = {region: color for region, color in anno_order_color.items() if region in distribution_percentage.columns}
## Plot the cumulative bar plot with custom colors
ax = distribution_percentage.plot(kind='bar', stacked=True, color=colors_for_species.values(), figsize=(10, 10), width=0.95)
## Customize the plot
plt.title('Dissection Distribution Across Types (Composition %)')
plt.xlabel(anno_level)
plt.ylabel('Composition Percentage (%)')
plt.legend(title='Dissection', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
# Save the plot
plt.savefig(os.path.join(figure_dir, f"cumulative_barplot_region_composition.pdf"))

###########################---------------------  (2-1) Cell counts per Group in each species  ---------------------  ###########################

## Define the colors for each subclass
bar_colors = anno.color_hex_group.tolist()

## Compute max count across all species and all categories
global_max = 0
# for species in adata_obs.organism.unique():
species_obs = adata_obs#[adata_obs['organism'] == species]
anno_table = species_obs[[anno_level, 'anatomical_region']]
anno_table[anno_level] = anno_table[anno_level].astype(str)
anno_table[anno_level] = pd.Categorical(anno_table[anno_level], categories=ct_order, ordered=True)
group_counts = anno_table[anno_level].value_counts().reindex(ct_order, fill_value=0)
species_max = group_counts.max()
if species_max > global_max:
    global_max = species_max

## Add 10% buffer for better visualization
global_max = global_max * 1.4

# for species in adata_obs.organism.unique():
## Species data
species_obs = adata_obs #[adata_obs['organism'] == species]
anno_table = species_obs[[anno_level, 'anatomical_region']]
## Anno table
anno_table[anno_level] = anno_table[anno_level].astype(str).copy()
anno_table[anno_level] = pd.Categorical(anno_table[anno_level], categories=ct_order, ordered=True)
## Counts
group_counts = anno_table[anno_level].value_counts()
group_counts = group_counts.reindex(ct_order, fill_value=0)
## Plot!
plt.figure(figsize=(10, 2))
## ax = group_counts.plot(kind='bar', color=bar_colors)
ax = plt.gca()
## Plot as points
x = range(len(ct_order))
y = group_counts.values
ax.scatter(x, y, color=bar_colors, zorder=10)
## Style
plt.yscale('log')
plt.axhline(y=100, color='grey', linestyle='--', linewidth=0.5, zorder=1)
plt.axhline(y=1000, color='grey', linestyle='--', linewidth=0.5, zorder=1)
plt.axhline(y=10000, color='grey', linestyle='--', linewidth=0.5, zorder=1)
ax.set_ylim(1, global_max)
plt.title('Number of Cells in Each Group')
plt.xlabel(anno_level)
plt.ylabel('Number of Cells (log scale)')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, f"cell_count_barplot.pdf"))

###########################---------------------  (2-2) Cluster counts per Group in each species  ---------------------  ###########################

# for species in adata_obs.organism.unique():
## Species data
species_obs = adata_obs#[adata_obs['organism'] == species]
anno_table = species_obs[[anno_level, 'anatomical_region', 'cluster_id']]
## Anno table
anno_table[anno_level] = anno_table[anno_level].astype(str).copy()
anno_table[anno_level] = pd.Categorical(anno_table[anno_level], categories=ct_order, ordered=True)
## Calculate the number of unique clusters per group
cluster_counts = anno_table.groupby(anno_level)['cluster_id'].nunique()
cluster_counts = cluster_counts.reindex(ct_order, fill_value=0)
## Plot the bar graph for the number of clusters in each group with log scale
plt.figure(figsize=(30, 10))
ax = cluster_counts.plot(kind='bar', color=bar_colors)
plt.yscale('log')
# Add gray dashed lines
plt.axhline(y=5, color='grey', linestyle='--', linewidth=1)
plt.axhline(y=10, color='grey', linestyle='--', linewidth=1)
ax.set_ylim(0, 100)  # Add a buffer of 10% to the max value
# Customize the plot
plt.title('Number of Clusters in Each Group')
plt.xlabel(anno_level)
plt.ylabel('Number of Clusters')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, f"cluster_count_barplot.pdf"))
