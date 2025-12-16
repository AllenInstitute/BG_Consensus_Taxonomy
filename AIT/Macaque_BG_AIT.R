library(scrattch.taxonomy)
library(anndata)
library(reticulate)

##
wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/s3/AIT"
setwd(wkdir)

## ------------------------------------------------------
adata = read_h5ad("Macaque_HMBA_basalganglia_AIT.h5ad");

##
isValid = checkTaxonomy(adata, print.messages=TRUE)

## Fixes
adata$obs$self_reported_ethnicity_ontology_term_id = adata$obs$self_reported_ethnicity
adata$obs$assay_ontology_term_id = "EFO:0030059"
adata$obs$organism_ontology_term_id = "NCBITaxon:9483"
adata$obs$anatomical_region = "Brain"
adata$obs$anatomical_region_ontology_term_id = "UBERON:0000955"
adata$obs$brain_region_ontology_term_id = "DHBA:10155"

adata$var$ensembl_id = NULL

adata$uns$mode = "standard"
adata$uns$filter = list()
adata$uns$filter["standard"] = 'False'

adata$uns$title = "Macaque_HMBA_basalganglia_consensus_AIT"
adata$uns$schema_version = as.character(packageVersion("scrattch.taxonomy"))

adata$uns$cluster_info = adata$obs[!duplicated(adata$obs$cluster_id),]

##
isValid = checkTaxonomy(adata, print.messages=TRUE)

## Finalize AIT
adata$write_h5ad("Macaque_HMBA_basalganglia_AIT_pre-print.h5ad")