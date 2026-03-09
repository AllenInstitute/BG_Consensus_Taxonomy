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
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/TAC3/"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/TAC3"

## Ensure cxgdir exists
os.makedirs(cxgdir, exist_ok=True)

## ----------------------------------------
adata = ad.read_h5ad(os.path.join(work_dir, "TAC3_data_merged_ortholog.h5ad"))

## Filter to class for both species
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="UMIs", batch_key="organism_sci", subset=False)
adata_hvg = adata[:,adata.var.highly_variable].copy()


genes_not_expressed = {}
for org in adata.obs["organism_sci"].unique():
    # Subset cells for this organism
    adata_sub = adata_hvg[adata.obs["organism_sci"] == org]
    # Count expression per gene across all cells
    gene_sums = np.array(adata_sub.X.sum(axis=0)).ravel()
    # Genes with zero counts
    zero_genes = adata_sub.var_names[gene_sums == 0]
    genes_not_expressed[org] = list(zero_genes)

# Convert to a dataframe for easier viewing
genes_not_expressed_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in genes_not_expressed.items()]))
genes_not_expressed_set = list(set(reduce(lambda x, y: x + y, genes_not_expressed.values())))

## ----------------------------------------
## Integrate across donors and species
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="organism_sci", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)

model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=128, 
                        n_latent=32, 
                        n_layers=2)
model.train(max_epochs=200)
 
## Save the latent space from scVI for downstream analysis
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata_hvg, use_rep = "X_scVI")
sc.tl.umap(adata_hvg, min_dist=0.3)
sc.tl.leiden(adata_hvg, flavor="igraph", n_iterations=2, resolution=0.3, key_added="leiden_X_scVI")

# ## Check for weird types in obs
# for col in adata_hvg.var.columns:
#     weird = (adata_hvg.var[[col]].applymap(type) != adata_hvg.var[[col]].iloc[0].apply(type)).any(axis=1)
#     if len(adata_hvg.var[weird]) > 0:
#         print(col)
#         adata_hvg.var[col] = adata_hvg.var[col].astype(str)
#     if col.endswith("ortholog_table_idx"):
#         adata_hvg.var[col] = adata_hvg.var[col].astype(str)
#     if col.endswith("ncbi_id"):
#         adata_hvg.var[col] = adata_hvg.var[col].astype(str)

# ## Check for weird types in obs
# for col in adata_hvg.obs.columns:
#     weird = (adata_hvg.obs[[col]].applymap(type) != adata_hvg.obs[[col]].iloc[0].apply(type)).any(axis=1)
#     if len(adata_hvg.obs[weird]) > 0:
#         print(col)
#         adata_hvg.obs[col] = adata_hvg.obs[col].astype(str)

# adata_hvg.obs_names_make_unique()
adata_hvg.write(cxgdir + "/TAC3_data_merged_ortholog_integrated.h5ad")

## ----------------------------------------
# adata_trim = ~adata.obs.leiden_X_scVI.isin(['15', '63', '66', '67', '69', '70', '71', '72', '73'])
# adata_trim = adata_trim & ~((adata.obs["AIT21.nbhd"] == "NN-IMN-GC") & adata.obs["macosko_fig1_region_label"].isin(["BS", "CB", "HPF", "HY", "Isocortex", "MY", "TH", "MB"]))
# adata_trim = adata_trim & ~(adata.obs["AIT21.subclass"] == "213 SCsg Gabrr2 Gaba")
# adata_trimmed = adata[adata_trim]
# adata_trimmed.write(work_dir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v3.h5ad")

## ----------------------------------------
# adata_trim = ~adata.obs["AIT21.nbhd"].isin(["Subpallium-GABA;HY-EA-Glut-GABA", "TH-EPI-Glut"])
# adata = adata[adata_trim]
# adata.write(work_dir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v6.h5ad")


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


