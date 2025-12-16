library(scrattch.taxonomy)
library(anndata)
library(reticulate)

##
wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/s3"
setwd(wkdir)

## ------------------------------------------------------
adata = read_h5ad("Macaque_basalganglia_pre-print_clean_pre-AIT_float32.h5ad");

## Ensure count matrix and annotations are in the same order.
taxonomy.counts = Matrix::t(adata$layers["UMIs"])
taxonomy.anno = adata$obs
taxonomy.anno$cluster = taxonomy.anno$cluster_id

## Compute binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster_id, 2000)

## Compute UMAP coordinates
umap.coords = adata$obsm$X_umap
rownames(umap.coords) = colnames(taxonomy.counts)

## Build Shiny taxonomy
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
                            meta.data = taxonomy.anno,
                            feature.set = binary.genes,
                            umap.coords = umap.coords,
                            title="HMBA_BG_consensus_macaque_AIT",
                            taxonomyDir="/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/s3",
                            hierarchy = c("Neighborhood", "Class", "Subclass", "Group", "cluster_id"),
                            subsample=1000,
                            dend=NULL)

## Add markers to dendrogram
# addDendrogramMarkers(AIT.anndata = AIT.anndata)

## 
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

## Provide hierarchy of the taxonomy
hierarchy = list("Neighborhood_label", "Class_label", "Subclass_label", "Group_label", "cluster_label")

## Build MapMyCells stats into AIT file for hierarchy mapping
AIT.anndata = addMapMyCells(AIT.anndata, hierarchy)
