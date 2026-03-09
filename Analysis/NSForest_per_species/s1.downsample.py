import anndata as ad
import numpy as np
import os

# ----------------------------
# Downsampling function
def downsample_adata(adata, group_key, max_cells=10000, random_state=0):
    np.random.seed(random_state)
    keep_idx = []

    for group, idx in adata.obs.groupby(group_key).groups.items():
        idx = np.array(idx)
        if len(idx) > max_cells:
            sampled = np.random.choice(idx, max_cells, replace=False)
        else:
            sampled = idx
        keep_idx.append(sampled)

    keep_idx = np.concatenate(keep_idx)
    return adata[keep_idx].copy()

# ----------------------------
# Configuration
species_list = ["Macaque", "Marmoset", "Human"]
anno_key = "cluster_id"
base_input_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/BICAN-releases/final"
base_output_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/BICAN-releases/final/downsample"

for species in species_list:
    print(f"Processing {species}...")

    # Input/output path
    input_file = os.path.join(base_input_dir, f"{species}_HMBA_basalganglia_AIT_pre-print.h5ad")
    output_file = os.path.join(base_output_dir, f"{species}_downsampled_{anno_key}_10k.h5ad")

    # Load and downsample
    adata = ad.read_h5ad(input_file)
    adata_ds = downsample_adata(adata, group_key=anno_key, max_cells=10000)

    # Save result
    adata_ds.write(output_file)
    print(f"Saved: {output_file}")

