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

## Load species and ancestral .h5ad
adata_files = {
    "species": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER.h5ad")),
    "ancestral": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER_ancestral.h5ad"))
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
phylo_df["branch_length_norm"] = phylo_df["branch_length"] / phylo_df["branch_length"].sum()


## --------------------------------------------
## Evo-distance metric on >= 90% alignments
## --------------------------------------------

## Setup alignment and evo distances (branch lengths)
alignment_matrix = merged_adata[:,phylo_df["node_name"]].layers[">=90%_aligned"].toarray()  # peaks × nodes
branch_lengths = phylo_df["branch_length_norm"].to_numpy()

## Compute evo distances
merged_adata.obs["evo_distance"] = (alignment_matrix * branch_lengths).sum(axis=1)

evo_dist_df = merged_adata.obs["evo_distance"]

## Plot distribution of evo distances
plt.figure(figsize=(6, 3))

sns.histplot(
    data=merged_adata.obs,
    x="evo_distance",
    hue="species",       # color by species
    bins=100,
    element="step",      # for clear overlay lines
    stat="density",
    common_norm=False,
    alpha=0.3
)

plt.xlabel("Evolutionary distance (branch-length weighted alignment)")
plt.ylabel("Density")
plt.title("Distribution of peak conservation across species")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "figures", "zoonomia_overlap_HALPER_ancestral_evodist.png"), dpi=900)

## Save
evo_dist_df.to_csv(os.path.join(analysis_dir, "analysis", "zoonomia_overlap_HALPER_ancestral_species_evodist.csv"), index=True)
