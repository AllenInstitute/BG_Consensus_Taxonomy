library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(anndata)
library(openxlsx)

##
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/figures/CRESted"
species_name = "human"

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
data = read.csv(file.path(figure_dir, "HAR_human_CRESted_contribution_scores_summary.csv"), row.names = 1)
# rownames(data) = data$simple_name
colnames(data) = gsub("_", " ", gsub("\\.", "-", colnames(data)))

## Load data
data_matrix = as.matrix(data)

## -------------------------------------------------
## Apply data normalization for plotting
## -------------------------------------------------

## min-max normalize along the columns
# data_matrix_minmax = apply(data_matrix, 1, function(x) {
#     (x - min(x)) / (max(x) - min(x))
# })

# ## Apply a which.max along the rows to get the max value for each cell type
# max_ct_idx = apply(data_matrix_minmax, 2, function(x) {
#     which.max(abs(x))
# })

# ## Find cell type with max expression per gene
# max_ct = rownames(data_matrix_minmax)[max_ct_idx]
# max_val = apply(data_matrix, 1, max)

# ## Combine into a data frame
# gene_table = data.frame(
#     gene = colnames(data_matrix_minmax),
#     max_cell_type = factor(max_ct, levels = cluster_meta$Group),
#     max_expr = max_val
# )

# ## Arrange genes by their max cell type
# gene_table_sorted = gene_table %>%
#     arrange(max_cell_type, desc(max_expr))

# ## ------------------------------------------------
# ## Draw basic heatmap for Figure 1B
# ## ------------------------------------------------

# ## From sorted gene table only keep the top 5 markers per cell type
# gene_table_sorted_filtered = gene_table_sorted %>%
#     group_by(max_cell_type) %>%
#     slice_head(n = 2) %>%
#     ungroup() %>%
#     as.data.frame()

data_matrix_plot = data_matrix[,-which(colnames(data_matrix) %in% c("BAM"))]

## ComplexHeatmap
heatmap = Heatmap(t(data_matrix_plot[names(sort(rowSums(abs(data_matrix_plot)), decreasing=T))[1:50],intersect(cluster_meta$Group, colnames(data_matrix_plot))]),
                cluster_columns = TRUE, 
                show_column_dend = FALSE,
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                col=circlize::colorRamp2(c(-0.5, -0.1, 0, 0.1, 0.5), c("blue", "white", "white", "white", "red")),
                # bottom_annotation = ha.tf_hit,
                # left_annotation = ha.celltype, 
                # row_split = celltype_sep,
                # column_split=col_sep,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title = NULL,
                border=TRUE,
                show_column_names = TRUE,
                show_row_names= TRUE)


png(paste0(file.path(figure_dir, "HAR_human_CRESted_celltypes_contributions.png")), width=10, height=5, units="in", res=600)
draw(heatmap)
dev.off()

