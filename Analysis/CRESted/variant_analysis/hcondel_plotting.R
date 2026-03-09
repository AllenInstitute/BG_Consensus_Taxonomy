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
data = read.csv(file.path(figure_dir, "hCONDEL_CRESted_predictions_difference.csv"))
rownames(data) = data$hCONDEL_ID
colnames(data) = gsub("\\.", "-", colnames(data))

## Load data
data_matrix = as.matrix(data[,-1])

## -------------------------------------------------
## Apply data normalization for plotting
## -------------------------------------------------

normalize_pos_neg <- function(x) {
  # Positive values
  pos <- x > 0
  if (any(pos)) {
    x[pos] <- x[pos] / (max(x[pos]))
  }
  # Negative values
  neg <- x < 0
  if (any(neg)) {
    x[neg] <- - x[neg] / (min(x[neg]))  # min(neg) is negative, so -x/min(neg) gives values in [-1,0]
  }
  return(x)
}

## min-max normalize along the columns
# data_matrix_minmax = matrix(0, nrow = nrow(adata$layers["mean"]), ncol = ncol(adata$layers["mean"]))
# data_matrix_minmax = apply(data_matrix, 1, function(x) {
#     (x - min(x)) / (max(x) - min(x))
# })

data_matrix_minmax <- t(apply(data_matrix, 1, normalize_pos_neg))

## Apply a which.max along the rows to get the max value for each cell type
# max_ct_idx = apply(adata$layers["mean_min-max"], 2, function(x) {
#     which.max(abs(x))
# })

# ## Find cell type with max expression per gene
# max_ct = adata$obs_names[max_ct_idx]
# max_val = apply(adata$layers["mean"], 2, max)

# ## Combine into a data frame
# gene_table = data.frame(
#     gene = adata$var_names,
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

data_matrix_minmax = data_matrix_minmax[-which(rowSums(data_matrix_minmax) == 0), ] 


## ComplexHeatmap
heatmap = Heatmap(t(data_matrix_minmax),
                cluster_columns = TRUE, 
                show_column_dend = FALSE,
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                col=colorRampPalette(c("red","lemonchiffon1","white","darkseagreen2", "blue"))(25),
                # bottom_annotation = ha.tf_hit,
                # left_annotation = ha.celltype, 
                # row_split = celltype_sep,
                # column_split=col_sep,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title = NULL,
                border=FALSE,
                show_column_names = FALSE,
                show_row_names= TRUE)


png(paste0(file.path(figure_dir, "hCONDEL_CRESted_human_model_differences.png")), width=10, height=5, units="in", res=600)
draw(heatmap)
dev.off()

