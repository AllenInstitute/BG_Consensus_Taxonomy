import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import glob 
from functools import reduce
import re
import os
from collections import Counter
import scvi

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"

cross_species_paths = {
    "HMBA:Human" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/Human_basalganglia_AIBS_BICAN.h5ad",
    "HMBA:Macaque" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/Macaque_HMBA_basalganglia_AIT_pre-print.h5ad",
    #"HMBA:Marmoset": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/marmoset/Marmoset_HMBA_basalganglia_AIT_pre-print.h5ad",
    #"Mouse_Broad": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/mouse/Broad_Mouse_BasalGanglia.h5ad"
}

## ----------------------------------------
## Load cross species data
cross_species = {}
for k in cross_species_paths.keys():
    adata = ad.read_h5ad(cross_species_paths[k])
    if k == "Mouse_Broad":
        adata.var_names = adata.var.gene_name
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        adata.layers["UMIs"] = adata.X.copy()  # Store the raw counts in a layer
        adata.obs["donor_id"] = "mouse"
        adata.obs["organism_sci"] = "Mus musculus"
        del adata.X
    elif k == "HMBA:Human":
        adata.layers["UMIs"] = adata.X.copy()  # Store the raw counts in a layer
        del adata.raw
        del adata.X
    else:
        adata.layers["UMIs"] = adata.raw.X.copy()  # Store the raw counts in a layer
        del adata.raw
        del adata.X
    cross_species[k] = adata

## ----------------------------------------
## Handle genes
ortholog_table = pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes2/all_gene_ids.csv', header=0)

species_conversion = {
    "Homo sapiens": "Homo_sapiens",
    "Macaca Mulatta": "Macaca_mulatta",
    "Callithrix jacchus": "Callithrix_jacchus",
    "Mus musculus": "Mus_musculus",
}

## Convert to orthologs, organism_ontology_term_id
for k in cross_species.keys():
    print(k)
    if k == "HMBA:Human":
        species = "Homo_sapiens"
        species_name = "Homo_sapiens"
    else:
        species = cross_species[k].obs['organism_sci'].mode()[0]
        species_name = species_conversion.get(species)
    ##
    if species_name is None:
        print(f"Species {species} not found in conversion dictionary for key: {k}")
        continue
    ##
    print("Identifying orthologs for: " + species,flush=True)
    ortholog_dict = dict(zip(
        ortholog_table.loc[ortholog_table['species'] == species_name, 'gene'], 
        ortholog_table.loc[ortholog_table['species'] == species_name, 'ortholog_symbol']
    ))
    ## Replace with orthologs
    cross_species[k].var[f'{species}_mgi_symbol'] = cross_species[k].var.index
    cross_species[k].var.index = cross_species[k].var.index.to_series().replace(ortholog_dict)
    cross_species[k].var_names = cross_species[k].var_names.astype('str')
    cross_species[k].var_names_make_unique()
    cross_species[k].var_names = cross_species[k].var_names.astype('category')
    ## Add tracking info
    tracker_dict = dict(zip(
        ortholog_table.loc[ortholog_table['species'] == species_name, 'gene'], 
        ortholog_table.loc[ortholog_table['species'] == species_name,:].index.tolist()
    ))
    cross_species[k].var[f'{species}_ortholog_table_idx'] = cross_species[k].var[f'{species}_mgi_symbol'].replace(tracker_dict)
    ## Add NCBI ID
    ncbi_dict = dict(zip(
        ortholog_table.loc[ortholog_table['species'] == species_name, 'gene'], 
        ortholog_table.loc[ortholog_table['species'] == species_name, 'ncbi_id']
    ))
    cross_species[k].var[f'{species}_ncbi_id'] = cross_species[k].var[f'{species}_mgi_symbol'].replace(ncbi_dict)

## Filter to common genes across species
common_genes = list(set.intersection(*map(set, [cross_species[key].var_names for key in cross_species.keys()])))
for dataset in cross_species.keys():
    print(dataset)
    cross_species[dataset] = cross_species[dataset][:,common_genes].copy()

## ----------------------------------------
## Merge and keep unique metadata from each object
adata = ad.concat(cross_species, join="outer", merge="unique")
del cross_species

## 
cols_to_drop = [col for col in adata.obs.columns if col.startswith('gex') or 
                                                    col.startswith('atac') or 
                                                    col.startswith('RNA') or 
                                                    col.startswith('ATAC')]
adata.obs.drop(columns=cols_to_drop, inplace=True)

## Filter to class for both species
# sc.pp.highly_variable_genes(adata_trimmed, flavor="seurat_v3", n_top_genes=3000, layer="UMIs", batch_key="organism_sci", subset=False)
# adata_hvg = adata_trimmed[:,adata_trimmed.var.highly_variable].copy()

marker_collection = []
for anno in ["neighborhood", "class", "subclass", "group"]:
    marker_set = pd.read_csv(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs/{anno}_tau_scores.csv", index_col=0)
    if anno != "group":
        marker_collection.append(marker_set.loc[marker_set.xspecies_min > 0.9].index.tolist())
    else:
        marker_collection.append(marker_set.loc[marker_set.xspecies_min > 0.75].index.tolist())

##
taxonomy_markers = list({item for sublist in marker_collection for item in sublist})

##
adata_hvg = adata[:, adata.var_names.isin(taxonomy_markers)].copy()

## ----------------------------------------
## Integrate across donors and species
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="organism", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
## Save model
model.save(dir_path=work_dir + "/pre-print-integration_human_base/")

## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
sc.tl.umap(adata, min_dist=0.3)
# sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=2, key_added="leiden_X_scVI")

##
# adata.X = adata.layers["UMIs"].copy()  # Use the UMIs layer for downstream analysis
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

## Check for weird types in obs
for col in adata.var.columns:
    weird = (adata.var[[col]].applymap(type) != adata.var[[col]].iloc[0].apply(type)).any(axis=1)
    if len(adata.var[weird]) > 0:
        print(col)
        adata.var[col] = adata.var[col].astype(str)
    if col.endswith("ortholog_table_idx"):
        adata.var[col] = adata.var[col].astype(str)
    if col.endswith("ncbi_id"):
        adata.var[col] = adata.var[col].astype(str)

## Save donor corrected latent space
# adata.write(cxgdir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v5.h5ad")

# ## Plot umap
# sc.pl.umap(adata, color=["organism_sci"], frameon=False, wspace=0.4, hspace=0.4, ncols=3, save="_HMBA_Human_Macaque_Marmoset_Mouse_BG_alignment.png")
# sc.pl.umap(adata, color=["leiden_X_scVI"], save="_HMBA_Human_Macaque_Marmoset_Mouse_BG_alignment_leiden.png")


adata = adata[:,['GATA3', 'PAX8']]
adata.obs = adata.obs.loc[:,["Group", "organism", "load_id", "donor_id"]]
adata.write(cxgdir + "/HMBA_Human_base_Macaque_BG_alignment.h5ad")

## ----------------------------------------
# adata_trim = ~adata.obs.leiden_X_scVI.isin(['15', '63', '66', '67', '69', '70', '71', '72', '73'])
# adata_trim = adata_trim & ~((adata.obs["AIT21.nbhd"] == "NN-IMN-GC") & adata.obs["macosko_fig1_region_label"].isin(["BS", "CB", "HPF", "HY", "Isocortex", "MY", "TH", "MB"]))
# adata_trim = adata_trim & ~(adata.obs["AIT21.subclass"] == "213 SCsg Gabrr2 Gaba")
# adata_trimmed = adata[adata_trim]
# adata_trimmed.write(cxgdir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v3.h5ad")
# adata_trimmed.write(work_dir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v3.h5ad")



## ----------------------------------------
# with h5py.File(cxgdir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v5.h5ad") as f:
#     adata_obs = ad.experimental.read_elem(f['obs'])
#     adata_obsm = ad.experimental.read_elem(f['obsm'])

# ## Add in mapping results from MapMyCells to primates
# ait_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/"
# label_df = pd.DataFrame()
# for species in ["Human", "Macaque", "Marmoset"]:
#     print(species)
#     if species in ["Human", "Macaque"]:
#         with h5py.File(f"{ait_dir}/{species}/BasalGanglia/{species}_basalganglia_AIBS_BICAN_MapMyCells.h5ad") as f:
#             species_obs = ad.experimental.read_elem(f['obs'])
#     elif species == "Marmoset":
#         with h5py.File(f"{ait_dir}/{species}/BasalGanglia/{species}_basalganglia_anno_latest_MapMyCells.h5ad") as f:
#             species_obs = ad.experimental.read_elem(f['obs'])
#         species_obs = species_obs.loc[:, ~species_obs.columns.str.contains(r'Flat|Siletti_|Siletti-')]
#     filtered_columns = [col for col in species_obs.columns if col.startswith(('Siletti', 'ABCmouse', 'HMBA_WB'))]
#     label_df = pd.concat([label_df, species_obs.loc[species_obs.index.isin(adata_obs.index), filtered_columns]], axis=0)
        
# ## Now merge adata_obs and label_df based on index
# merged_obs = adata_obs.merge(label_df, left_index=True, right_index=True, how='left')














## ----------------------------------------
adata = ad.read_h5ad(cxgdir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v5.h5ad")
adata_trimmed = ~adata.obs["AIT21.nbhd"]










# ##
# import scanpy as sc
# import scipy
# adata = sc.read_h5ad(datadir + "macaque_mouse_species_integrated_SpC.h5ad")

# ## Get centroids for each subclass in scVI integrated space
# homology_df = {
#     "macaque": pd.DataFrame(index=adata.obs.subclass.unique(), columns=["centroid"]),
#     "mouse": pd.DataFrame(index=adata.obs.subclass.unique(), columns=["centroid"]),
#     "homology": pd.DataFrame(index=adata.obs.subclass.unique(), columns=adata.obs.subclass.unique())
# }

# ## Compute centroids per species
# for subclass in adata.obs.subclass.unique():
#     homology_df["macaque"]["centroid"][subclass] = np.median(adata[(adata.obs.subclass == subclass) & (adata.obs.species == "Macaca nemestrina")].obsm["X_umap_species_integrated"], axis=0)
#     homology_df["mouse"]["centroid"][subclass] = np.median(adata[(adata.obs.subclass == subclass) & (adata.obs.species == "Mouse")].obsm["X_umap_species_integrated"], axis=0)

# ## Build homology matrix
# for subclass in adata.obs.subclass.unique():
#     for subclass_query in adata.obs.subclass.unique():
#         if np.isnan(homology_df["macaque"]["centroid"][subclass]).any() | np.isnan(homology_df["mouse"]["centroid"][subclass_query]).any():
#             homology_df["homology"][subclass][subclass_query] = 0
#         else:
#             homology_df["homology"][subclass][subclass_query] = scipy.spatial.distance.euclidean(homology_df["macaque"]["centroid"][subclass], homology_df["mouse"]["centroid"][subclass_query])

# # ##
# # import pickle 
# # with open('homology_mouse_macaque_SpC.pkl', 'wb') as f:
# #     pickle.dump(homology_df, f)

# # ##
# # with open('homology_mouse_macaque_SpC.pkl', 'rb') as f:
# #     homology_df = pickle.load(f)

# ##
# homology_df["homology"].to_csv("homology_mouse_macaque_SpC.csv")

##
# library(ComplexHeatmap)
# library(circlize)

# ##
# mat = read.csv("homology_mouse_macaque_SpC.csv")
# mat.plot = as.matrix(mat[,-1]); rownames(mat.plot) = colnames(mat.plot) = mat[,1]
# mat.plot = mat.plot[rowSums(mat.plot) != 0, colSums(mat.plot) != 0]

# ##
# # mat.plot = mat.plot[rownames(mat.plot)[apply(mat.plot, 2, which.min)],]
# # mat.plot = mat.plot[,colnames(mat.plot)[apply(mat.plot, 1, which.min)]]

# ##
# heatmap = Heatmap(mat.plot,
#                   cluster_columns = T, 
#                   cluster_rows = T,
#                   col=colorRamp2(breaks=c(0, 1, 2), colors=c("black", "grey", "white")),
#                   # top_annotation = ha.top,
#                   # left_annotation = ha.row, 
#                   #row_split=gene.anno$set_anno,
#                   #row_km=4,
#                   #column_split=col_sep,
#                   row_names_gp = gpar(fontsize = 8),
#                   column_names_gp = gpar(fontsize = 8),
#                   border=T,
#                   show_column_names=T,
#                   show_row_names=T)

# pdf(paste0("~/mouse_macaque_SpC_homology.pdf"), width=7, height=7)
# draw(heatmap)
# dev.off()
