## singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.6
library(scrattch.taxonomy)
library(scrattch.mapping)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

## ------------------------------------------------------
## Load in the marmoset .h5ad after the alignment + preprocessing + QC steps
adata_marmoset = read_h5ad("./HMBA_BG_Marmoset_raw.h5ad")

## ------------------------------------------------------
## Map against Human 
AIT.anndata = loadTaxonomy(taxonomyDir="/PATH/TO/TAXONOMY/Human/",
                            anndata_file="HMBA_Human_BG_082024_AIT.h5ad")

## Map! Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata ,
                                query.data=Matrix::t(adata_marmoset$X),
                                label.cols=c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label"), ## Which obs in AIT.anndata contain annotations to map. E.g. "class", "subclass", etc.
                                corr.map=TRUE,
                                tree.map=FALSE,
                                hierarchical.map=TRUE,
                                seurat.map=FALSE)

## Extract mapping results from S4 mappingClass
mapping.results = list()
mapping.results[["Human"]] = getMappingResults(mapping.anno, scores = TRUE)
write.csv(mapping.results[["Human"]], file="./mapping_results_Human.csv")

##
colnames(mapping.results[["Human"]]) = paste0("HMBA_Human_Consensus-", colnames(mapping.results[["Human"]]))
adata_marmoset$obs = cbind(adata_marmoset$obs, mapping.results[["Human"]])

## ------------------------------------------------------
## Map against Macaque
AIT.anndata = loadTaxonomy(taxonomyDir="/PATH/TO/TAXONOMY/Macaque/",
                            anndata_file="HMBA_Macaque_BG_082024_AIT.h5ad")

## Map! Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata ,
                                query.data=Matrix::t(adata_marmoset$X),
                                label.cols=c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label"), ## Which obs in AIT.anndata contain annotations to map. E.g. "class", "subclass", etc.
                                corr.map=TRUE,
                                tree.map=FALSE,
                                hierarchical.map=TRUE,
                                seurat.map=FALSE)

## Extract mapping results from S4 mappingClass
mapping.results[["Macaque"]] = getMappingResults(mapping.anno, scores = FALSE)
write.csv(mapping.results[["Macaque"]], file="./mapping_results_Macaque.csv")

##
colnames(mapping.results[["Macaque"]]) = paste0("HMBA_Macaque_Consensus-", colnames(mapping.results[["Macaque"]]))
adata_marmoset$obs = cbind(adata_marmoset$obs, mapping.results[["Macaque"]])

## ------------------------------------------------------
## Save
adata_marmoset$write_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Regional_Taxonomies/BasalGanglia/Marmoset/HMBA_BG_Marmoset_mapped.h5ad")