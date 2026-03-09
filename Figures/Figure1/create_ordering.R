## Load scrattch.taxonomy, use scrattch docker!
## docker://njjai/scrattch_mapping:1.0.0
library(scrattch.taxonomy)
library(reticulate)
library(anndata)
library(dendextend)
cell_type_mapper <- import("cell_type_mapper")

##
# taxonomy = read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/xspecies/Consensus_HMBA_basalganglia_AIT_pre-print.h5ad", backed="r")

# ## Group level marker genes
# marker_files = list.files("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs/", pattern="*_tau_scores.csv", full.names=TRUE)
# marker_collection = list()
# for(file in marker_files) {
#     markers = read.csv(file, row.names=1)
#     if(grepl("group", file)){
#       marker_collection[[file]] = markers %>% filter(markers$xspecies_min > 0.75) %>% rownames()
#     } else {
#       marker_collection[[file]] = markers %>% filter(markers$xspecies_min > 0.9) %>% rownames()
#     }
# }
# marker_set = Reduce(union, marker_collection)

# ## Filter taxonomy for marker genes
# taxonomy = taxonomy[,taxonomy$var_names %in% marker_set]

# ##
# cluster = taxonomy$obs$Group; 
# cluster = gsub(" ", "_", cluster) # Replace spaces with underscores
# names(cluster) = taxonomy$obs_names

# ##
# expr.data = Matrix::t(taxonomy$layers["UMIs"])
# expr.data = as(expr.data, "dgCMatrix")
# expr.data = scrattch.bigcat::logCPM(as.matrix(expr.data))

# ##
# cluster_stats = scrattch.bigcat::get_cl_medians(expr.data, cluster)
# cluster_stats = cluster_stats[-which(rowSums(cluster_stats) == 0),]

# ##
# write.csv(cluster_stats, "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/dendrogram/cluster_stats.csv")

cluster_stats = read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/dendrogram/cluster_stats.csv", row.names=1)

##
celltype_order = c(
  "STRd D1 Matrix MSN",
  "STRd D1 Striosome MSN",
  "STRd D2 Matrix MSN",
  "STRd D2 Striosome MSN",
  "STRv D1 MSN",
  "STRv D2 MSN",
  "STRd D2 StrioMat Hybrid MSN",
  "STR D1D2 Hybrid MSN",
  "STRv D1 NUDAP MSN",
  "OT D1 ICj",
  "GPe MEIS2-SOX6 GABA",
  "AMY-SLEA-BNST GABA",
  "AMY-SLEA-BNST D1 GABA",
  "SN SEMA5A GABA",
  "SN EBF2 GABA",
  "OB FRMD7 GABA",
  "OB Dopa-GABA",
  "GPi Shell",
  "GPi Core",
  "SN GATA3-PVALB GABA",
  "GPe SOX6-CTXND1 GABA",
  "ZI-HTH GABA",
  "SN-VTR-HTH GATA3-TCF7L2 GABA",
  "SN GATA3-PAX8 GABA",
  "GPe-NDB-SI LHX6-LHX8-GBX1 GABA",
  "VIP GABA",
  "LAMP5-CXCL14 GABA",
  "STR SST-CHODL GABA",
  "STR SST-RSPO2 GABA",
  "STR LYPD6-RSPO2 GABA",
  "STR SST-ADARB2 GABA",
  "LAMP5-LHX6 GABA",
  "STR FS PTHLH-PVALB GABA",
  "STR TAC3-PLPP4 GABA",
  "STR-BF TAC3-PLPP4-LHX8 GABA",
  "GPin-BF Cholinergic GABA",
  "STR Cholinergic GABA",
  "STRd Cholinergic GABA",
  "SN-VTR GAD2 Dopa",
  "SN SOX6 Dopa",
  "SN-VTR CALB1 Dopa",
  "BF SKOR1 Glut",
  "VTR-HTH Glut",
  "STH PVALB-PITX2 Glut",
  "OPC",
  "COP",
  "ImOligo",
  "Oligo PLEKHG1",
  "Oligo OPALIN",
  "ImAstro",
  "Astrocyte",
  "Ependymal",
  "Microglia",
  "BAM",
  "Monocyte",
  "B cells",
  "T cells",
  "VLMC",
  "SMC",
  "Pericyte",
  "Endo"
)

##
celltype_order = gsub(" ", "_", celltype_order)
l.rank = setNames(1:length(celltype_order), celltype_order)

##
colnames(cluster_stats) = gsub("\\.", "-", colnames(cluster_stats))

##
dend.result = build_dend(
    cl.dat  = cluster_stats,
    l.rank  = l.rank, 
    nboot   = 2,
    ncores  = 1)

## Plot
dend = dend.result$dend

## Customize the dendrogram using dendextend
dend_plot <- dend %>%
  # set("branches_k_color", k=3) %>%   # Color branches based on 3 clusters
  set("labels_cex", 0.8) %>%           # Adjust label size
  set("branches_lwd", 1)               # Adjust branch thickness
dend_plot <- set(dend_plot, "branches_col", "black")

## Plot the dendrogram
pdf("/home/nelson.johansen/dendrogram.pdf", width=15, height=7)
par(mar = c(15, 4, 4, 2))
plot(dend_plot, main = "Enhanced Dendrogram with dendextend")
dev.off()

## Extract the order of labels from the dendrogram
dend_order <- labels(dend_plot)
dend_order <- gsub("_", " ", dend_order) 

dend_order = data.frame(
    group = dend_order,
    order = seq_along(dend_order)
)

## Add order to metadata table
anno_table = read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation - consensus_anno_pre-print.csv")

##
anno_table$display_order_group = anno_table[match(anno_table$Group, dend_order$group), "display_order_group"]

write.csv(anno_table, 
    "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/group_order_anno_pre-print.csv",
    row.names=FALSE)