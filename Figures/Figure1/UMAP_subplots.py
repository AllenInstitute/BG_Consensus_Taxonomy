import os
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

##
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/figures"

##
adata = ad.read_h5ad(os.path.join(data_dir, "Consensus_HMBA_basalganglia_AIT_pre-print.h5ad"), backed="r")

## 
metadata = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
colors = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="value_set_colors")
colors_map = dict(zip(colors["label"], colors["color_hex_triplet"]))

group_color_map = dict(zip(metadata["Group"], metadata["color_hex_group"]))

## --------------------------------------------------
## UMAP colored by Group
## --------------------------------------------------

group_colors = adata.obs["Group"].map(group_color_map)

# Plot
plt.figure(figsize=(20, 20))
scatter = plt.scatter(adata.obsm["X_umap"][:,0], 
                    adata.obsm["X_umap"][:,1], 
                    c=group_colors, 
                    s=0.05, 
                    alpha=0.8)

plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP colored by cell type")
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, "figure1", "UMAP_by_group_HMBA_BG_consensus.png"), dpi=400)

## --------------------------------------------------
## UMAP colored by region (Human + Macaque) and grey (Marmoset)
## --------------------------------------------------
adata_region = adata[adata.obs["organism"].isin(["Human", "Macaque"])] ## Marmoset has no anatomical region info
adata_region_marmoset = adata[adata.obs["organism"].isin(["Marmoset"])] ## Marmoset has no anatomical region info

## Coarsen anatomical regions
region_map = {
    "NACc": "NAC",
    "NACs": "NAC",
    "GPeR": "GPe",
    "GPeC": "GPe",
    "SN-VTA": "SN",
}
region_labels = adata_region.obs["anatomical_region"].replace(region_map)
region_colors = region_labels.map(colors_map)

# Plot
plt.figure(figsize=(20, 20))
plt.scatter(adata_region_marmoset.obsm["X_umap"][:,0], 
                adata_region_marmoset.obsm["X_umap"][:,1], 
                c=["#D3D3D3"] * adata_region_marmoset.n_obs, ## Grey for Marmoset "Brain"
                s=0.05, 
                alpha=0.8)
plt.scatter(adata_region.obsm["X_umap"][:,0], 
                adata_region.obsm["X_umap"][:,1], 
                c=region_colors, 
                s=0.05, 
                alpha=0.8)
##
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP colored by cell type")
plt.tight_layout()
plt.savefig(os.path.join(figure_dir, "figure1", "UMAP_by_region_HMBA_BG_consensus.png"), dpi=400)

