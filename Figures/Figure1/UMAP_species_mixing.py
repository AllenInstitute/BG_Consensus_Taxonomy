import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sciduck as sd
from scipy.stats import entropy
from tqdm import tqdm
from anndata import AnnData
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import matplotlib as mpl


##
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/figures"

## Compute species entropy per cell
def cell_entropy(adata: AnnData,
                    annotation_columns: list,
                    nearest_neighbors: int = 15, 
                    dim: str = "X_scVI"
                ) -> AnnData:
    """
    Compute entropy (mixing) of annotations within a cells local neighborhood. 

    Args:
        adata: AnnData object with `annotations`
        annotation_columns: Cell level annotations.
        nearest_neighbors: Number of nearest neighbors.
        dim: Dimensionality reduction
        
    Returns:
        Returns the updated AnnData object.
    """
    def neighbors_from_adata(adata, n_neighbors=15):
        if "distances" not in adata.obsp:
            print("Computing neighbors...")
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=dim)
        dists = adata.obsp["distances"]
        indices = np.array([row.indices[:n_neighbors] for row in dists.tocsr()])
        return indices
    ## Build nearest neighbor tree for fast lookup
    nearest_ind = neighbors_from_adata(adata, nearest_neighbors)
    ##
    for anno in annotation_columns:
        print("Computing entropy on: " + anno)
        adata.obs[anno + "_entropy"] = -1 ## Initialize with a value outside range of entropy.
        for cell in tqdm(range(0, adata.shape[0])):
            nearest_neighbors = nearest_ind[cell,:]
            adata.obs.loc[adata.obs.index[cell], anno + "_entropy"] = entropy(adata.obs.loc[adata.obs.index[nearest_neighbors],anno].value_counts()/len(nearest_neighbors))
    return adata

##
adata = ad.read_h5ad(os.path.join(data_dir, "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"), backed="r")

##
adata = cell_entropy(adata, annotation_columns=["organism"], nearest_neighbors=100, dim="X_scVI")

##
perm = np.random.permutation(adata.n_obs)

## Plot
plt.figure(figsize=(20, 20))
scatter = plt.scatter(adata.obsm["X_umap"][perm,0], 
                adata.obsm["X_umap"][perm,1], 
                cmap="magma",
                c=adata.obs.organism_entropy[perm], ## Grey for Marmoset "Brain"
                s=0.05, 
                alpha=0.8)
##
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP colored by cell type")
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, "figure1", "UMAP_species_mixing_HMBA_BG_consensus_legend.png"), dpi=400)

## Draw the legend separately
cmap = plt.get_cmap("magma")
norm = mpl.colors.Normalize(vmin=adata.obs.organism_entropy.min(),
                            vmax=adata.obs.organism_entropy.max())
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 

fig, ax = plt.subplots(figsize=(2, 5))  # narrow figure for vertical colorbar
cbar = fig.colorbar(sm, cax=ax)          # pass the axes explicitly
cbar.set_label("Organism Entropy")

# Save as separate PNG
plt.savefig(os.path.join(figure_dir, "figure1", "UMAP_entropy_legend.png"), dpi=400, bbox_inches="tight")
plt.close()