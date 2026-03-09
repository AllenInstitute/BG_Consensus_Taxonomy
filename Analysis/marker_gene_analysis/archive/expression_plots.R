library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(anndata)
library(openxlsx)

base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package"
analysis_dir = file.path(base_dir, "analysis", "expressologs")

## Save 
figure_dir = file.path(analysis_dir, "figures")
if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE)
}

## AnnoTable
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet="consensus_anno_pre-print")
cluster_meta = cluster_meta %>%
    distinct(Group, .keep_all = TRUE)
cluster_meta = cluster_meta[cluster_meta$display_order_group,]

## Plot for all data objects
data = "AIBS_BG_consensus_AIT-pre-print_conserved_markers.h5ad" # "AIBS_BG_consensus_AIT-pre-print_tf_markers.h5ad" 
print(paste0("Processing: ", gsub(".h5ad", "", data)))

## Load data
adata = read_h5ad(file.path(analysis_dir, data))
adata = adata[cluster_meta$Group]

## -------------------------------------------------
## Apply data normalization for plotting
## -------------------------------------------------

## min-max normalize along the columns
adata$layers["mean_min-max"] = matrix(0, nrow = nrow(adata$layers["mean"]), ncol = ncol(adata$layers["mean"]))
adata$layers["mean_min-max"] = apply(adata$layers["mean"], 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
})

## Apply a which.max along the rows to get the max value for each cell type
max_ct_idx = apply(adata$layers["mean_min-max"], 2, function(x) {
    which.max(x)
})

## Find cell type with max expression per gene
max_ct = adata$obs_names[max_ct_idx]
max_val = apply(adata$layers["mean"], 2, max)

## Combine into a data frame
gene_table = data.frame(
    gene = adata$var_names,
    max_cell_type = factor(max_ct, levels = cluster_meta$Group),
    max_expr = max_val
)

## Arrange genes by their max cell type
gene_table_sorted = gene_table %>%
    arrange(max_cell_type, desc(max_expr))

## ------------------------------------------------
## Draw basic heatmap for Figure 1B
## ------------------------------------------------

## From sorted gene table only keep the top 5 markers per cell type
gene_table_sorted_filtered = gene_table_sorted %>%
    group_by(max_cell_type) %>%
    slice_head(n = 2) %>%
    ungroup() %>%
    as.data.frame()

##
adata_plot = adata[,gene_table_sorted_filtered$gene]

## ComplexHeatmap
heatmap = Heatmap(adata_plot$layers["mean_min-max"],
                cluster_columns = FALSE, 
                show_column_dend = FALSE,
                cluster_rows = FALSE,
                col=colorRampPalette(c("white","darkseagreen2","darkseagreen2", "blue"))(25),
                # bottom_annotation = ha.tf_hit,
                # left_annotation = ha.celltype, 
                # row_split = celltype_sep,
                # column_split=col_sep,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title = NULL,
                border=FALSE,
                show_column_names = FALSE,
                show_row_names= FALSE)


png(paste0(file.path(figure_dir, paste0(gsub(".h5ad", "", data), "_marker_plot_group_figure1b.png"))), width=5, height=20, units="in", res=600)
draw(heatmap)
dev.off()


## ------------------------------------------------
## Draw heatmap split by Class for visualization of complete marker sets
## ------------------------------------------------

##
adata_plot = adata[,gene_table_sorted$gene]

## Create cell type blocks and annotation bar
class_meta = cluster_meta %>%
    distinct(Class, .keep_all = TRUE)
class_colors = class_meta$color_hex_class
names(class_colors) = class_meta$Class

ha.celltype = rowAnnotation(
    celltype = factor(cluster_meta$Class, levels = class_meta$Class),
    col = list(celltype = class_colors),
    show_annotation_name = FALSE
)

## Heatmap blocks
celltype_sep = factor(cluster_meta$Class, levels = class_meta$Class)

## ComplexHeatmap
heatmap = Heatmap(adata_plot$layers["mean_min-max"],
                cluster_columns = FALSE, 
                show_column_dend = FALSE,
                cluster_rows = FALSE,
                col=colorRampPalette(c("white","darkseagreen2","darkseagreen2", "blue"))(25),
                # bottom_annotation = ha.tf_hit,
                left_annotation = ha.celltype, 
                row_split = celltype_sep,
                # column_split=col_sep,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title = NULL,
                border=TRUE,
                show_column_names = TRUE,
                show_row_names= TRUE)

## Save 
figure_dir = file.path(analysis_dir, "figures")
if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE)
}

png(paste0(file.path(figure_dir, paste0(gsub(".h5ad", "", data), "_marker_plot_group.png"))), width=10, height=5, units="in", res=600)
draw(heatmap)
dev.off()
