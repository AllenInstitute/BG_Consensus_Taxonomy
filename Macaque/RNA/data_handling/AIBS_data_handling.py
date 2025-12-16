import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import os

##
species = "Macaque"

##
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/bicore/OCS_multiome_datasets/macaque/BG_paper_Macaque"

##
Q19_Q21 = ad.read_h5ad(os.path.join(work_dir, "Q19-Q21", "Q19-Q21_gex.h5ad"))
Q23 = ad.read_h5ad(os.path.join(work_dir, "Q23", "Q23_gex.h5ad"))

##
adata = ad.concat([Q19_Q21, Q23], merge="same", uns_merge="same")

##
adata.write(os.path.join(work_dir, f"{species}_basalganglia_AIBS.h5ad"))

# ##
species = "Macaque"
data_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
# adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_NeMO.h5ad"))

##
adata = ad.read_h5ad("Macaque_basalganglia_NeMO.h5ad")
mean_expression = np.array(adata.X.mean(axis=0)).flatten()
gene_names = adata.var_names 
mean_expression_df = pd.DataFrame({'Gene': gene_names, 
                                        'Mean Expression': mean_expression})
mean_expression_df.to_csv(os.path.join(work_dir, f"{species}_basalganglia_NeMO_mean_expression.csv"), index=False)

adata = ad.read_h5ad(os.path.join(work_dir, "Macaque_basalganglia_AIBS.h5ad"))
mean_expression = np.array(adata.X.mean(axis=0)).flatten()
gene_names = adata.var_names 
mean_expression_df = pd.DataFrame({'Gene': gene_names, 
                                        'Mean Expression': mean_expression})
mean_expression_df.to_csv(os.path.join(work_dir, f"{species}_basalganglia_AIBS_mean_expression.csv"), index=False)
