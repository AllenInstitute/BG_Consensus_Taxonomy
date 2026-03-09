import anndata as ad
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/analysis/"

## human cactus matrices
custom_adata = ad.read_h5ad(os.path.join(work_dir, "human_zoonomia_overlap.h5ad"))
halper_adata = ad.read_h5ad(os.path.join(work_dir, "human_zoonomia_overlap_HALPER.h5ad"))

# Compute means for each variable
custom_means = np.asarray(custom_adata.X.mean(axis=0)).flatten()
halper_means = np.asarray(halper_adata.X.mean(axis=0)).flatten()

##
plt.figure(figsize=(8, 6))
sns.kdeplot(custom_means, label="Custom", color="blue", fill=True, alpha=0.4)
sns.kdeplot(halper_means, label="HALPER", color="orange", fill=True, alpha=0.4)

plt.title("Density Curves of Mean Values")
plt.xlabel("Mean Value")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig("/home/nelson.johansen/mean_density_comparison.png", dpi=300)