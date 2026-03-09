import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

##
##
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## AnnoTable
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Read data
species_h5ads = {
    "human" : os.path.join(data_dir, "human", "crested_adata", "human_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "macaque" : os.path.join(data_dir, "macaque", "crested_adata", "macaque_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "marmoset": os.path.join(data_dir, "marmoset", "crested_adata", "marmoset_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
}

## Load epi-conserved-marker orthologs
epi_conserved_markers = pd.read_csv(os.path.join(analysis_dir, "atac", "conservation", "epi_conserved_marker_peaks.csv"))

species_hits = {}
for species, h5ad_path in species_h5ads.items():
    print(f"Processing {species}...")
    species_adata = ad.read_h5ad(h5ad_path)
    species_adata = species_adata[:, species_adata.var_names.isin(epi_conserved_markers[f'{species}_ID'])]
    species_adata.var["gini_scores"] = epi_conserved_markers.set_index(f"{species}_ID").loc[species_adata.var_names, f"{species}_gini"].values
    print(f"Number of peaks after filtering to level 3 orthologs: {species_adata.n_vars}")
    ## Identify most accessible cell type for each peak
    max_idx = np.argmax(species_adata.X, axis=0)
    max_cell_types = species_adata.obs_names[max_idx]
    species_adata.var["max_cell_type"] = max_cell_types.tolist()
    ## Gather top 100 peaks per cell type based on specificity (Gini)
    top_peaks = []
    celltype_record = pd.DataFrame(columns=["region", "cell_type"])
    for cell_type in species_adata.obs_names.unique():
        cell_type_peaks = species_adata.var[species_adata.var["max_cell_type"] == cell_type]
        top_cell_type_peaks = cell_type_peaks.nlargest(10000, f"gini_scores")
        top_peaks.append(top_cell_type_peaks.index.tolist())
        ##
        celltype_df = pd.DataFrame({
            "region": top_cell_type_peaks.index,
            "cell_type": [cell_type] * top_cell_type_peaks.shape[0]
        })
        celltype_record = pd.concat([celltype_record, celltype_df], axis=0)
    ##
    top_peaks = [peak for sublist in top_peaks for peak in sublist]  # Flatten list
    top_peaks = list(set(top_peaks))  # Unique peaks
    ## Record 
    species_hits[species] = celltype_record
    ## Subset for heatmap
    heatmap_adata = species_adata.copy()#[:,top_peaks].copy()
    heatmap_adata.obs["cell_type"] = heatmap_adata.obs_names.astype("category")
    ## Set nan to 0
    heatmap_adata.X = np.nan_to_num(heatmap_adata.X, nan=0.0)
    ## Apply min-max normalization to the columns for better visualization
    heatmap_adata.X = (heatmap_adata.X - heatmap_adata.X.min(axis=0)) / ((heatmap_adata.X.max(axis=0) - heatmap_adata.X.min(axis=0)) + 1e-9)
    ##
    del heatmap_adata.uns
    del heatmap_adata.obsm
    del heatmap_adata.layers
    ## Normalize for visualization
    # sc.pp.normalize_total(heatmap_adata, target_sum=1e4)
    # sc.pp.log1p(heatmap_adata)
    # sc.pp.scale(heatmap_adata)
    ## Arrange by max_cell_type var and Group order
    group_order = cluster_meta.Group
    ## If any Groups are missing add a obs of nans
    missing_groups = [grp for grp in group_order if grp not in heatmap_adata.obs_names]
    for missing in missing_groups:
        missing_adata = ad.AnnData(
            X=np.full((1, heatmap_adata.n_vars), np.nan),
            obs=pd.DataFrame(index=[missing]),
            var=heatmap_adata.var.copy()
        )
        missing_adata.obs["cell_type"] = missing
        heatmap_adata = ad.concat([heatmap_adata, missing_adata], axis=0, join="outer", merge="same")
    ## Ensure correct categorical ordering for obs
    heatmap_adata = heatmap_adata[heatmap_adata.obs_names.isin(group_order),:]
    heatmap_adata.obs["cell_type"] = heatmap_adata.obs["cell_type"].astype(
        pd.CategoricalDtype(categories=group_order, ordered=True)
    )
    ## Ensure correct categorical ordering by arranging var max_cell_type following group_order
    heatmap_adata = heatmap_adata[:, heatmap_adata.var_names[heatmap_adata.var["max_cell_type"].isin(group_order)]]
    heatmap_adata.var["max_cell_type"] = heatmap_adata.var["max_cell_type"].astype(
        pd.CategoricalDtype(categories=group_order, ordered=True)
    )
    heatmap_adata = heatmap_adata[:, heatmap_adata.var.sort_values("max_cell_type").index]
    if species == "human":
        peak_order = heatmap_adata.var_names
    else:
        heatmap_adata = heatmap_adata[:, peak_order.intersection(heatmap_adata.var_names)]
    ## Color
    cmap = copy.copy(cm.get_cmap("viridis"))
    cmap.set_bad(color="grey")  # NaNs will be grey
    ## Plot heatmap
    sc.pl.heatmap(
        heatmap_adata,
        var_names=heatmap_adata.var_names,
        groupby="cell_type", # Or another category for annotation
        show_gene_labels=False,
        dendrogram=False,
        figsize=(12,10),
        vmin=0,
        vmax=1,
        cmap=cmap,
        show=False
    )
    plt.savefig(os.path.join(analysis_dir,"atac", "conservation", "figures", f"{species}_level3_heatmap.png"), dpi=700)

# ## Save level 3 peaks
# level3_human = species_hits["human"]
# level3_human.to_csv(os.path.join(out_dir, "level3_ortholog_peaks_per_cell_type_human.csv"), index=False)

## Merge dataframes
all_top_peaks = pd.concat(species_hits, axis=0, names=['species']).reset_index(level='species').reset_index(drop=True)
all_top_peaks.to_csv(os.path.join(analysis_dir,"atac", "conservation", "ortholog_peaks_per_cell_type_per_species.csv"), index=False)

