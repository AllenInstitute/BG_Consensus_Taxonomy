## --------------------------------------------------
## Plot phylogenetic tree of species in cactus HAL
## --------------------------------------------------
import os
import pandas as pd
import anndata as ad
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

##
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Load in peak x group boolean matrix
for species in ["marmoset"]:
    print(f"Processing {species}...")
    peak_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/{species}/ATAC"
    peak_group = pd.read_csv(os.path.join(peak_dir, "Group_by_peaks.csv"), index_col=0)
    peak_group.columns = [col.replace("_", " ") for col in peak_group.columns]
    ## Compute jaccard co-accessibility on peak_group
    cell_types = cluster_meta.Group.tolist()
    n_ct = len(cell_types)
    ##
    jaccard = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)
    cond_prob = pd.DataFrame(index=cell_types, columns=cell_types, dtype=float)
    for i in peak_group.columns:
        vec_i = peak_group[i].values.astype(bool)
        n_i = vec_i.sum()
        for j in peak_group.columns:
            vec_j = peak_group[j].values.astype(bool)
            ##
            both = (vec_i & vec_j).sum()
            either = (vec_i | vec_j).sum()
            ## Compute Jaccard and conditional probabilities
            jaccard.loc[i, j] = both / either if either > 0 else np.nan
            cond_prob.loc[i, j] = both / n_i if n_i > 0 else np.nan
    ## --------------------------------------------
    ## Color
    cmap = copy.copy(cm.get_cmap("Blues"))
    cmap.set_bad(color="lightgrey")  # NaNs will be grey
    ## Plot heatmap with seaborn
    plt.figure(figsize=(20, 20))
    sns.heatmap(
        jaccard,
        cmap=cmap,
        vmin=0,
        vmax=1,
        square=True,
        linewidths=0,
        linecolor="none",
        cbar=False,  # remove colorbar
        mask=jaccard.isna()  # ensures NaNs are grey
    )
    plt.title(f"{species} - Cell-type Jaccard co-accessibility", fontsize=14)
    plt.xlabel("Cell type")
    plt.ylabel("Cell type")
    plt.tight_layout()
    ##
    plt.savefig(
        os.path.join(
            analysis_dir, "atac", "conservation", "figures", f"{species}_coaccessibility_jaccard_heatmap.png"
        ),
        dpi=900
    )
    plt.close()