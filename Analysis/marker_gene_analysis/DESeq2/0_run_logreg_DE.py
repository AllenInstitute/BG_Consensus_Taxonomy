import anndata as ad
import scanpy as sc
import xarray as xr
import pandas as pd
import numpy as np
import os

## Helper
def _rgg_df_all_groups(species_adata, groupby, key="rank_genes_groups"):
    try:
        # works in newer Scanpy
        return sc.get.rank_genes_groups_df(species_adata, None, key=key)
    except TypeError:
        # older Scanpy: must request one group at a time
        frames = []
        for g in list(species_adata.obs[groupby].cat.categories):
            try:
                df_g = sc.get.rank_genes_groups_df(species_adata, g, key=key)
                frames.append(df_g)
            except KeyError:
                continue
        if not frames:
            return pd.DataFrame(columns=["names","group","scores"])
        return pd.concat(frames, ignore_index=True)


## Data location
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/de_genes_scanpy/logreg"

## -------------------------------------------
## Identify species files, could also be a list from glob
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "Human_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "Macaque_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "Marmoset_HMBA_basalganglia_AIT_pre-print_orthologs.h5ad"),
}

## Parameters
top_k = 25
min_cells_per_group = 5
groupby = "Group"
species_key = "organism_name"
layer = None
use_raw = False

## Collect group levels and genes
group_levels = set()
genes = set()

species_score_mats = {}
for sp in species_h5ads.keys():
    species_adata = ad.read_h5ad(species_h5ads[sp])
    ## Ensure categorical
    if not pd.api.types.is_categorical_dtype(species_adata.obs[groupby]):
        species_adata.obs[groupby] = species_adata.obs[groupby].astype("category")
    group_levels.update(species_adata.obs[groupby].cat.categories)
    genes.update(species_adata.var_names.to_list())
    ## Filter groups with too few cells
    counts = species_adata.obs[groupby].value_counts()
    keep_groups = counts[counts >= min_cells_per_group].index.tolist()
    if len(keep_groups) < 2:
        species_score_mats[sp] = pd.DataFrame(np.nan, index=genes, columns=group_levels)
        continue
    species_adata = species_adata[species_adata.obs[groupby].isin(keep_groups)].copy()
    species_adata.obs[groupby] = species_adata.obs[groupby].astype("category")
    ## Run logreg
    sc.tl.rank_genes_groups(
        species_adata,
        groupby=groupby,
        method="logreg",
        n_genes=species_adata.n_vars,
        layer=layer,
        use_raw=use_raw
    )
    ## Extract results
    df = _rgg_df_all_groups(species_adata, groupby, key="rank_genes_groups")
    df.to_csv(os.path.join(analysis_dir, f"{sp}_logreg_all_groups.csv"), index=False)
    # if df.empty:
    #     scmat = pd.DataFrame(np.nan, index=genes, columns=group_levels)
    # else:
    #     scmat = (df.pivot(index="names", columns="group", values="scores")
    #                 .reindex(index=genes, columns=group_levels))
    # species_score_mats[sp] = scmat

# data = np.stack([species_score_mats[sp].to_numpy().T for sp in species_levels], axis=0)
# da = xr.DataArray(
#     data,
#     coords={"species": species_levels, "group": group_levels, "gene": genes},
#     dims=("species", "group", "gene"),
#     name="logreg_score"
# )
