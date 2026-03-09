library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(anndata)
library(openxlsx)

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/figures/"

## AnnoTable
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet="consensus_anno_pre-print")
cluster_meta = cluster_meta %>%
    distinct(Group, .keep_all = TRUE)
cluster_meta = cluster_meta[cluster_meta$display_order_group,]
cluster_meta$groupby = cluster_meta$Group

## Files to plot
files = c("AIBS_BG_consensus_taxonomy_harmonized_AIT-pre-print_conserved_markers_10925.h5ad", 
          "AIBS_BG_consensus_taxonomy_harmonized_AIT-pre-print_tf_markers_10925.h5ad")

## AIBS_SpC_consensus_taxonomy_harmonized_AIT-pre-print_tf_markers_091725.h5ad
for(file in files){ 

    print(paste0("Processing: ", gsub(".h5ad", "", file)))

    ## Load data
    adata = read_h5ad(file.path(analysis_dir, "marker_gene_analysis", file)) # AIBS_SpC_consensus_taxonomy_harmonized_AIT-pre-print_conserved_markers.h5ad
    adata = adata[cluster_meta$Group]

    ## Create cell type blocks and annotation bar
    class_meta = cluster_meta %>%
        distinct(Neighborhood, .keep_all = TRUE)
    class_colors = class_meta$color_hex_neighborhood
    names(class_colors) = class_meta$Neighborhood

    ## Bring in TF analysis from CRESted
    # tf_crested_hits = read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/CRESted/tf_crested_conserved.csv")
    # tf_crested_hits = tf_crested_hits %>%
    #     filter(TF %in% adata$var_names) 

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
    max_ct <- adata$obs_names[max_ct_idx]
    max_val <- apply(adata$layers["mean"], 2, max)

    # Combine into a data frame
    gene_table <- data.frame(
        gene = adata$var_names,
        #tau = adata$var[""],
        max_cell_type = factor(max_ct, levels = cluster_meta$Group),
        max_expr = max_val
    )

    # Arrange genes by their max cell type
    gene_table_sorted <- gene_table %>%
        arrange(max_cell_type, desc(max_expr))

    ## From sorted gene table only keep the top markers per cell type
    gene_table_sorted_filtered = gene_table_sorted %>%
        group_by(max_cell_type) %>%
        slice_head(n = 10) %>%
        ungroup() %>%
        as.data.frame()

    ## 
    gene_table_sorted_filtered$Neighborhood = cluster_meta$Neighborhood[match(gene_table_sorted_filtered$max_cell_type, cluster_meta$Group)]
    gene_table_sorted_filtered$Class = cluster_meta$Class[match(gene_table_sorted_filtered$max_cell_type, cluster_meta$Group)]

    ##
    adata_plot = adata[,gene_table_sorted_filtered$gene]

    ## Create 
    # tf_hits_human = adata_plot$var_names %in% (tf_crested_hits %>% filter(species == "Human") %>% pull(TF))
    # tf_hits_macaque = adata_plot$var_names %in% (tf_crested_hits %>% filter(species == "Macaque") %>% pull(TF))
    # tf_hits_mouse = adata_plot$var_names %in% (tf_crested_hits %>% filter(species == "Mouse") %>% pull(TF))
    # tf_hits_combined = adata_plot$var_names %in% (tf_crested_hits %>% filter(species == "Combined") %>% pull(TF))

    # ha.tf_hit <- HeatmapAnnotation(
    #     Human     = anno_simple(as.character(tf_hits_human), col = c("TRUE"="black","FALSE"="white"), height = unit(1.25, "mm")),
    #     Macaque   = anno_simple(as.character(tf_hits_macaque), col = c("TRUE"="black","FALSE"="white"), height = unit(1.25, "mm")),
    #     Mouse     = anno_simple(as.character(tf_hits_mouse), col = c("TRUE"="black","FALSE"="white"), height = unit(1.25, "mm")),
    #     Conserved = anno_simple(as.character(tf_hits_combined), col = c("TRUE"="black","FALSE"="white"), height = unit(1.25, "mm")),
    #     show_legend = FALSE,
    #     show_annotation_name = TRUE,
    #     annotation_name_gp = gpar(fontsize = 6),
    #     annotation_name_side = "left",
    #     border = TRUE
    # )

    adata_plot$obs$Neighborhood = cluster_meta$Neighborhood[match(adata_plot$obs_names, cluster_meta$Group)]
    adata_plot$obs$Class = cluster_meta$Class[match(adata_plot$obs_names, cluster_meta$Group)]

    ## Heatmap blocks
    adata_plot_neuron = adata_plot[adata_plot$obs$Neighborhood != "Nonneuron",]
    celltype_sep = factor(adata_plot_neuron$obs$Neighborhood, levels = class_meta$Neighborhood)

    ha.celltype <- HeatmapAnnotation(
        celltype = factor(adata_plot_neuron$obs$Neighborhood, levels = class_meta$Neighborhood),
        col = list(celltype = class_colors),
        show_annotation_name = FALSE, 
        show_legend = FALSE
    )

    ## Genes
    gene_table_sorted_filtered %>% filter(Neighborhood != "Nonneuron") %>% pull(gene) -> gene_set

    ## ComplexHeatmap Neurons
    heatmap = Heatmap(t(adata_plot_neuron[,gene_set]$layers["mean_min-max"]),
                    cluster_columns = FALSE, 
                    show_column_dend = FALSE,
                    cluster_rows = FALSE,
                    col=colorRampPalette(c("white","darkseagreen2","darkseagreen2", "blue"))(25),
                    # bottom_annotation = ha.tf_hit,
                    top_annotation = ha.celltype, 
                    column_split = celltype_sep,
                    # column_split=col_sep,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_title = NULL,
                    column_title = NULL,
                    border=TRUE,
                    show_heatmap_legend = FALSE,
                    show_column_names = TRUE,
                    show_row_names= TRUE)

    png(file.path(figure_dir, paste0(gsub(".h5ad", "", file), "_group_neurons.png")), width=10, height=8, units="in", res=600)
    draw(heatmap)
    dev.off()

    ## ComplexHeatmap Non-Neurons
    adata_plot_nonneuron = adata_plot[adata_plot$obs$Neighborhood == "Nonneuron",]
    celltype_sep = factor(adata_plot_nonneuron$obs$Neighborhood, levels = class_meta$Neighborhood)

    ha.celltype <- HeatmapAnnotation(
        celltype = factor(adata_plot_nonneuron$obs$Neighborhood, levels = class_meta$Neighborhood),
        col = list(celltype = class_colors),
        show_annotation_name = FALSE,
        show_legend = FALSE
    )

    ## Genes
    gene_table_sorted_filtered %>% filter(Neighborhood == "Nonneuron") %>% pull(gene) -> gene_set

    heatmap = Heatmap(t(adata_plot_nonneuron[,gene_set]$layers["mean_min-max"]),
                    cluster_columns = FALSE, 
                    show_column_dend = FALSE,
                    cluster_rows = FALSE,
                    col=colorRampPalette(c("white","darkseagreen2","darkseagreen2", "blue"))(25),
                    # bottom_annotation = ha.tf_hit,
                    top_annotation = ha.celltype, 
                    column_split = celltype_sep,
                    # column_split=col_sep,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_title = NULL,
                    column_title = NULL,
                    border=TRUE,
                    show_heatmap_legend = FALSE,
                    show_column_names = TRUE,
                    show_row_names= TRUE)

    png(file.path(figure_dir, paste0(gsub(".h5ad", "", file), "_group_non-neurons.png")), width=10, height=8, units="in", res=600)
    draw(heatmap)
    dev.off()

}