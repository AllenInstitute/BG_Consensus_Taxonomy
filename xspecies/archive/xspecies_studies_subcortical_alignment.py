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
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/xspecies"

cross_species_paths = {
    "HMBA:Human" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/Human_basalganglia_HMBA_AIT19-5_anno_latest.h5ad",
    "HMBA:Macaque" : "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia/Macaque_basalganglia_HMBA_AIT11-9_anno_latest.h5ad",
    "HMBA:Marmoset": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Marmoset/SubCortical/rna_subcortex_marm.h5ad"
}

# /allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Mouse/supplemental_data/Mouse_basalganglia_anno_latest.h5ad

cross_species = {k:sc.read_h5ad(v) for k,v in cross_species_paths.items()}

## Set layers before removing genes
cross_species["HMBA:Human"].X = cross_species["HMBA:Human"].raw.X.copy()
cross_species["HMBA:Macaque"].X = cross_species["HMBA:Macaque"].raw.X.copy()
cross_species["HMBA:Marmoset"].X = cross_species["HMBA:Marmoset"].layers["UMIs"].copy()

del cross_species["HMBA:Human"].raw
del cross_species["HMBA:Macaque"].raw
del cross_species["HMBA:Marmoset"].raw
del cross_species["HMBA:Marmoset"].layers["UMIs"]

## ----------------------------------------
## Handle genes
ortholog_table = pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes2/all_gene_ids.csv', header=0)

species_conversion = {
    "Human": "Homo_sapiens",
    "Macaque": "Macaca_mulatta",
    "Marmoset": "Callithrix_jacchus",
    "Mouse": "Mus_musculus",
}

cross_species["HMBA:Marmoset"].obs["organism"] = "Marmoset"

## Convert to orthologs, organism_ontology_term_id
for k in cross_species.keys():
    species = cross_species[k].obs['organism'].mode()[0]
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
    ##
    cross_species[k].var['mgi_symbol'] = cross_species[k].var.index
    cross_species[k].var.index = cross_species[k].var.index.to_series().replace(ortholog_dict)
    cross_species[k].var_names_make_unique()

cross_species["HMBA:Human"].layers["UMIs"] = cross_species["HMBA:Human"].X.copy()
cross_species["HMBA:Macaque"].layers["UMIs"] = cross_species["HMBA:Macaque"].X.copy()
cross_species["HMBA:Marmoset"].layers["UMIs"] = cross_species["HMBA:Marmoset"].X.copy()

## Filter to common genes across species
common_genes = list(set.intersection(*map(set, [cross_species[key].var_names for key in cross_species.keys()])))
for dataset in cross_species.keys():
    cross_species[dataset] = cross_species[dataset][:,common_genes].copy()

## -------------------
## Adjust metadata

cross_species["HMBA:Human"].obs.rename(columns={"cluster_id": "Human_cluster"}, inplace=True)
cross_species["HMBA:Macaque"].obs.rename(columns={"cluster_id": "Macaque_cluster"}, inplace=True)
cross_species["HMBA:Marmoset"].obs.rename(columns={"Cluster_tc_marm": "Marmoset_cluster"}, inplace=True)

def extract_terms(elements):
    return [e.split("_")[2] if len(e.split("_")) > 2 else None for e in elements]

cross_species["HMBA:Marmoset"].obs["donor_id"] = "unknown"
cross_species["HMBA:Marmoset"].obs["donor_id"] = extract_terms(cross_species["HMBA:Marmoset"].obs.index.tolist())
cross_species["HMBA:Marmoset"].obs["donor_id"] = cross_species["HMBA:Marmoset"].obs["donor_id"].astype('category')

## ----------------------------------------
## Merge and keep unique metadata from each object
adata = ad.concat(cross_species, join="outer")
del cross_species

## Filter to class for both species
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=6000, layer="UMIs", batch_key="organism", subset=False)
adata_hvg = adata[:,adata.var.highly_variable].copy()

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
# model.save(dir_path=work_dir + "/multi_study_integration_models/")

## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI_xspecies"] = model.get_latent_representation()
 
## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI_xspecies")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=5, key_added="leiden_scVI")
sc.tl.umap(adata, min_dist=0.3)

##---
adata.obsm["X_umap_species_integrated"] = adata.obsm["X_umap"].copy()

##
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

##
# to_keep = ['donor_id', 'organism', "roi",
#             "Human_cluster", "Macaque_cluster", "Marmoset_cluster", 
#             "Group", "Subclass", "Class", "Neighborhood",
#             "Group_2024_marm", "Group_scanvi_2024_marm"]
# adata.obs = adata.obs.loc[:,to_keep]

## Save donor corrected latent space
adata.write(cxgdir + "/HMBA_BG_Human_Macaque_Marmoset_alignment.h5ad")
adata.write(work_dir + "/HMBA_BG_Human_Macaque_Marmoset_alignment.h5ad")

# ## Write out the leiden_scVI table for quality control
# anno_table = sd.anno.build_annotation_table(adata, 
#                                     group_by="leiden_scVI", 
#                                     categorical_annotations=["donor_id"],
#                                     numeric_annotations=["doublet_score", "total_genes", "total_counts"], 
#                                     min_percent=0.05, 
#                                     annotation_alerts={"donor_name": 0.90})
# anno_table.to_csv(os.path.join(work_dir, f"alignment_leiden_consensus_anno_table.csv"))

# ###############################################
# ## Map labels across species / studies
# ###############################################
# ## Load in macaque spc reference
# mac_ref = ad.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/Macaque/data/NHP_SpC_AIT12_taxonomy.h5ad")

# ## Bring in macaque labels to aligned object
# adata.obs.loc[mac_ref.obs.index, "Class"] = mac_ref.obs["Class"]
# adata.obs.loc[mac_ref.obs.index, "Subclass"] = mac_ref.obs["Subclass"]
# adata.obs.loc[mac_ref.obs.index, "Group"] = mac_ref.obs["Group"]

# ## Build KNN classifier of scVI aligned space to transfer labels 
# from sklearn.neighbors import KNeighborsClassifier
# neigh=KNeighborsClassifier(n_neighbors=15)

# ## Transfer Class labels
# neigh.fit(adata[mac_ref.obs.index].obsm['X_scVI'], adata[mac_ref.obs.index].obs['Class'])
# annotated_labels = neigh.predict(adata[~adata.obs.index.isin(mac_ref.obs.index)].obsm['X_scVI'])
# adata.obs.loc[~adata.obs.index.isin(mac_ref.obs.index), "Class"] = annotated_labels

# ## Transfer Subclass labels
# neigh.fit(adata[mac_ref.obs.index].obsm['X_scVI'], adata[mac_ref.obs.index].obs['Subclass'])
# annotated_labels = neigh.predict(adata[~adata.obs.index.isin(mac_ref.obs.index)].obsm['X_scVI'])
# adata.obs.loc[~adata.obs.index.isin(mac_ref.obs.index), "Subclass"] = annotated_labels

# ## Transfer Group labels
# neigh.fit(adata[mac_ref.obs.index].obsm['X_scVI'], adata[mac_ref.obs.index].obs['Group'])
# annotated_labels = neigh.predict(adata[~adata.obs.index.isin(mac_ref.obs.index)].obsm['X_scVI'])
# adata.obs.loc[~adata.obs.index.isin(mac_ref.obs.index), "Group"] = annotated_labels

# ## Save
# adata.write(cxgdir + "multispecies_integrated_SpC_aligned_labels_realigned.h5ad")

























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
