import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
from scipy.sparse import csr_matrix

##
wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/"
os.chdir(wkdir)

## ----------------------------------------
## Load data
cross_species = {
    "macaque" : sc.read_h5ad("./HMBA_Macaque_BG_082024_AIT.h5ad"),
    "human": sc.read_h5ad("./HMBA_Human_BG_082024_AIT.h5ad"),
    "marmoset" : sc.read_h5ad("./marm_hmba_h5ad/marmoset_HMBA_BG_merged.h5ad"),
}

## ----------------------------------------
## Handle genes
ortholog_table = pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes/all_gene_ids.csv', header=0)

## Loop through species and update gene names to orthologs
for species, species_name in zip(['macaque', 'human', 'marmoset'], 
                                 ['Macaca_nemestrina', "Homo_sapiens", "Callithrix_jacchus"]):
    ##
    print("Identifying orthologs for: " + species)
    ortholog_dict = dict(zip(ortholog_table.loc[ortholog_table['species']==species_name,'gene'], 
                             ortholog_table.loc[ortholog_table['species']==species_name,'ortholog_symbol']))
    cross_species[species].var['mgi_symbol'] = cross_species[species].var.index
    cross_species[species].var.index = cross_species[species].var.index.to_series().replace(ortholog_dict)
    cross_species[species].var_names_make_unique()

## ----------------------------------------
## Overlap and filter genes
common_labels = list(set.intersection(*map(set, [cross_species[key].var_names for key in cross_species.keys()])))
for dataset in cross_species.keys():
    cross_species[dataset] = cross_species[dataset][:,common_labels].copy()

## Add UMI layer, ensure sparse format for scVI
cross_species["macaque"].X = csr_matrix(cross_species["macaque"].raw[:,common_labels].X.copy())
cross_species["human"].X = csr_matrix(cross_species["human"].raw[:,common_labels].X.copy())
cross_species["marmoset"].X = cross_species["marmoset"].raw[:,common_labels].X.copy()

## Add consistent Donor metadata column
cross_species["marmoset"].obs["Donor_label"] = cross_species["marmoset"].obs.donor_name

## Just ignore the raw data, subsetting is annoying
del cross_species["macaque"].raw
del cross_species["human"].raw
del cross_species["marmoset"].raw

## ----------------------------------------
## Merge and keep unique metadata from each object
adata = cross_species["macaque"].concatenate(cross_species["human"], cross_species["marmoset"], join="inner", batch_key="Study")
adata.layers["UMIs"] = adata.X.copy()

## ----------------------------------------
## Run scVI
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000, layer="UMIs", batch_key="Study", subset=False)

##
adata_hvg = adata[:,adata.var.highly_variable].copy()

## Integrate across donors and species
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="Study", ## We will model the gene dispersion per species
    categorical_covariate_keys=["Donor_label"],
)
model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=300)

## Save model
model.save(dir_path="HMBA_integration_BG")

## Extract integrated space
adata.obsm["X_scVI"] = model.get_latent_representation()

## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep="X_scVI")
adata.obsm["X_umap_integrated"] = sc.tl.umap(adata, min_dist=0.3, copy=True)
adata.obsm["X_umap_integrated"] = adata.obsm["X_umap_integrated"].obsm["X_umap"]

## ----------------------------------------
## Update metadata that is inconsistent (NA + Str in same field, etc.)
for col in adata.obs.columns:
    weird = (adata.obs[[col]].applymap(type) != adata.obs[[col]].iloc[0].apply(type)).any(axis=1)
    if len(adata.obs[weird]) > 0:
        print(col)
        adata.obs[col] = adata.obs[col].astype(str)

## ----------------------------------------
## Save
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/cross-species/basal_ganglia_HMBA_integrated.h5ad")
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalNuclei/basal_ganglia_HMBA_integrated.h5ad")

####################################################################################################
## Basic scanpy processing
####################################################################################################

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_scVI")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

####################################################################################################
## Prune Marmoset data to match Macaque and Human HMBA Basal Ganglia consensus taxonomies
####################################################################################################

## Compute Study entropy
for anno in ["Study"]:
    adata.obs[anno + "_entropy"] = -1 ## Initialize with a value outside range of entropy.
    for cluster in np.unique(adata.obs.leiden):
        adata.obs.loc[adata.obs.leiden == cluster, anno + "_entropy"] = scipy.stats.entropy(adata.obs.loc[adata.obs.leiden == cluster, anno].value_counts()/sum(adata.obs.leiden == cluster))

## Remove all clusters from only one study (Marmoset, in this case) using threshold of 0.1 on study entropy
adata = adata[adata.obs["Study_entropy"] > 0.1].copy()

## Save
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/cross-species/basal_ganglia_HMBA_integrated_pruned.h5ad")
adata.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalNuclei/basal_ganglia_HMBA_integrated_pruned.h5ad")

####################################################################################################
## Extract cleaned Marmoset data
####################################################################################################

## Remove study from index
adata.obs.index = [re.sub("-2$", "", s) for s in  adata.obs.index] ## remove study suffix, annoying

## Remove low quality cluster
adata = adata[~(adata.obs["leiden"] == '1')] ## Low quality cluster based on qc-metrics

## Load back in the original marmoset data and subset to common index between aligned and original data
adata_marmoset = sc.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Marmoset/raw_data/marm_hmba_h5ad/marmoset_HMBA_BG_merged.h5ad")
adata_marmoset = adata_marmoset[adata_marmoset.obs.index.isin(adata.obs.index)]

##
adata_marmoset.write("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Marmoset/HMBA_BG_Marmoset_raw.h5ad")


