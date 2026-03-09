import scanpy as sc
import numpy as np

import milopy
import milopy.core as milo
import milopy.utils

## Load in backed
adata = sc.read("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalNuclei/siletti_human_nhp_species_integrated_complete.h5ad", backed="r")

##
adata.obs.load_name = adata.obs.load_name.astype(str)
adata.obs.loc[adata.obs.study == '1', "load_name"] = adata.obs.loc[adata.obs.study == '1', "sample_id"].astype(str)

## Build KNN graph
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_scVI")

## Assign cells to neighbourhoods
milo.make_nhoods(adata)

## Count cells from each sample in each nhood
milo.count_nhoods(adata, sample_col="load_name")

## Test for differential abundance between conditions
milo.DA_nhoods(adata, design="~ study")

# ##
# milopy.utils.annotate_nhoods(adata, anno_col='cluster')

## Check results
milo_results = adata.uns["nhood_adata"].obs
milo_results

##
milopy.utils.build_nhood_graph(adata,
                                basis = "X_umap_species_integrated")

##
with rc_context({"figure.figsize": (9, 8), "figure.dpi": (600)}):
    plot_nhood_graph_2(adata, 
                            alpha=0.01, ## SpatialFDR level (1%) 
                            min_size=1, ## Size of smallest dot
                            save = "milo_nhood_graph.png",
                        )


##
def plot_nhood_graph_2(
    adata: AnnData,
    alpha: float = 0.1,
    min_logFC: float = 0,
    min_size: int = 10,
    plot_edges: bool = False,
    title: str = "DA log-Fold Change",
    **kwargs
):
    nhood_adata = adata.uns["nhood_adata"].copy()
    if "Nhood_size" not in nhood_adata.obs.columns:
        raise KeyError(
            'Cannot find "Nhood_size" column in adata.uns["nhood_adata"].obs -- \
                please run milopy.utils.build_nhood_graph(adata)'
        )
    nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"]
                        > alpha, "graph_color"] = np.nan
    nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
    nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"]
                        < min_logFC, "graph_color"] = np.nan
    # Plotting order - extreme logFC on top
    nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(),
                        "abs_logFC"] = np.nan
    ordered = nhood_adata.obs.sort_values(
        'abs_logFC', na_position='first').index
    nhood_adata = nhood_adata[ordered]
    vmax = 2
    vmin = -2
    sc.pl.embedding(nhood_adata, "X_milo_graph",
                    color="graph_color", cmap="RdBu_r",
                    size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    # edge_width =
                    sort_order=False,
                    frameon=False,
                    vmax=vmax, vmin=vmin,
                    title=title,
                    **kwargs
                    )