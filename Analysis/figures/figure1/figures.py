import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## Load HMBA consensus data
adata = ad.read_h5ad(os.path.join(work_dir, "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"))

## Neurotransmitter plots
adata.obs["neurotransmitter"] = ""
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("Glut"), "Glutamatergic", adata.obs["neurotransmitter"])
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("GABA"), "GABAergic", adata.obs["neurotransmitter"])
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("MSN"), "GABAergic", adata.obs["neurotransmitter"])
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("Dopa"), "Dopaminergic", adata.obs["neurotransmitter"])
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("Chol"), "Cholinergic", adata.obs["neurotransmitter"])
adata.obs["neurotransmitter"] = np.where(adata.obs["Subclass"].str.contains("GABA-Glut"), "GABAergic/Glutamatergic", adata.obs["neurotransmitter"])

## Plot neurotransmitter UMAP
plt.figure(figsize=(12,12))
neurotransmitter_palette = {
    "Glutamatergic": "#1f77b4",
    "GABAergic": "#ff7f0e",
    "Dopaminergic": "#2ca02c",
    "Cholinergic": "#d62728",
    "GABAergic/Glutamatergic": "#9467bd",
    "": "lightgrey"
}
plt.scatter(adata.obsm["X_umap"][:,0], 
            adata.obsm["X_umap"][:,1],
            c=adata.obs["neurotransmitter"].map(neurotransmitter_palette).values.astype(str),
            s=0.1,
            alpha=0.8)
plt.title(f"UMAP of HMBA BG Neurons - Neurotransmitter", fontsize=16)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "figure1", f"hmba_neuron_neurotransmitter_umap.png"), dpi=900)

## Count the number of cells in each Group per species
group_counts = adata.obs.groupby(["organism", "Group"]).size().unstack(fill_value=0)

## normalize per Group
group_props = group_counts.div(group_counts.sum(axis=0), axis=1)

## order the Groups by cluster_meta
ordered_groups = cluster_meta['Group'].tolist()
group_props = group_props[ordered_groups]

## Plot as a 3 column (species) by Group heatmap
plt.figure(figsize=(10,20))
## Color as white to black
sns.heatmap(group_props.T, 
            cmap="Greys", 
            cbar_kws={'label': 'Proportion of cells'}, 
            linewidths=0.5, 
            linecolor='gray')
plt.title("Proportion of cells per Group by Species", fontsize=16)
plt.xlabel("Species")
plt.ylabel("Group")
plt.savefig(os.path.join(analysis_dir, "figure1", f"hmba_neuron_group_species_heatmap.png"), dpi=900)
plt.savefig(os.path.join(analysis_dir, "figure1", f"hmba_neuron_group_species_heatmap.pdf"), dpi=900)

## Count the total number of cells per Group across species
total_group_counts = adata.obs['Group'].value_counts().reindex(ordered_groups).fillna(0)    

## Log 10 transform
total_group_counts_log = np.log10(total_group_counts + 1)

## Plot heatmap of total counts
plt.figure(figsize=(5,20))
sns.heatmap(total_group_counts_log.to_frame(), 
            cmap="Reds", 
            cbar_kws={'label': 'Log10(Number of cells + 1)'}, 
            linewidths=0.5, 
            linecolor='gray')
plt.title("Total number of cells per Group (log10 transformed)", fontsize=16)
plt.xlabel("Group")
plt.ylabel("")
plt.savefig(os.path.join(analysis_dir, "figure1", f"hmba_neuron_total_group_counts_heatmap.png"), dpi=900)
plt.savefig(os.path.join(analysis_dir, "figure1", f"hmba_neuron_total_group_counts_heatmap.pdf"), dpi=900)