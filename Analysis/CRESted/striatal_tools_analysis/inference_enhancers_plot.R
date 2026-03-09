library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(anndata)
library(openxlsx)
library(viridisLite)

## -------------------------
##
## -------------------------
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = read.xlsx(file.path(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet="consensus_anno_pre-print")
cluster_meta = cluster_meta %>% distinct(Group, .keep_all=TRUE)
cluster_meta$Group = gsub(" ", "_", cluster_meta$Group)

## -------------------------
## Load data
## -------------------------
adata = read_h5ad(file.path(analysis_dir, "crested", "data", "striatal_tools_enhancer_embeddings.h5ad"))

## Create color map
color_map = data.frame("group" = cluster_meta$Group, "color" = cluster_meta$color_hex_group)
color_map = color_map %>% filter(group %in% adata$var_names)

## -------------------------
## Plot Heatmap
## -------------------------

## Normalize data matrix
norm_matrix = apply(adata$X, 1, function(x) {
    (x - min(x)) / (max(x) - min(x))
})
adata$layers["min-max"] = t(norm_matrix)

## Gather prediction values for max cell type per enhancer
max_ct_idx = apply(adata$layers["min-max"], 1, function(x) {
    which.max(x)
})
max_pred_value = adata$X[cbind(1:nrow(adata), max_ct_idx)]
max_pred_value = log(max_pred_value, base=2)

ha.max.pred <- rowAnnotation(
    magnitude = max_pred_value,
    col = list(magnitude = colorRamp2(c(min(max_pred_value), max(max_pred_value)), c("white", "blue"))),
    show_annotation_name = FALSE, 
    show_legend = FALSE
)

## Colors
breaks <- seq(0, 1, length.out = 100)
col_fun <- colorRamp2(breaks, viridis(length(breaks)))

## Enhancer groupings
enh_groups = factor(adata$obs[["Cell type"]], levels=rev(c("pan-MSN", "D1 MSN", "D2 MSN", "Cholinergic", "Pvalb-Pthlh", "Sst-Chodl")))

## ComplexHeatmap Neurons
heatmap = Heatmap(adata[,color_map$group]$layers["min-max"],
                cluster_columns = FALSE, 
                show_column_dend = FALSE,
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                col=col_fun,
                right_annotation = ha.max.pred,
                # bottom_annotation = ha.tf_hit,
                #top_annotation = ha.celltype, 
                #column_split = celltype_sep,
                row_split=enh_groups,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title = NULL,
                column_title = NULL,
                border=TRUE,
                show_heatmap_legend = FALSE,
                show_column_names = TRUE,
                show_row_names= TRUE)

pdf(file.path(analysis_dir, "crested", "figures", "striatal_tools_enhancers.pdf"), width=10, height=8)
draw(heatmap)
dev.off()