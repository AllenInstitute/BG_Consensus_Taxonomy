library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(anndata)
library(openxlsx)

base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/"
analysis_dir = file.path(base_dir, "Analysis", "species_specific_markers")

## Save 
figure_dir = file.path(analysis_dir, "figures")
if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE)
}

## AnnoTable
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta = cluster_meta %>%
    distinct(Group, .keep_all = TRUE)
cluster_meta = cluster_meta[cluster_meta$group_hierarchy_diagram_order,]

color_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet="value_set_colors")
color_meta = color_meta %>%
    filter(field == "species_genus")

## Plot for all data objects
adata_set = list()
for(species in c("Human", "Macaque", "Mouse")) {

    ##
    print(paste0("Processing species: ", species))

    ## Gather data info
    data = paste0(tolower(species),"_exclusive_tau_markers.h5ad")
    print(paste0("Processing: ", gsub(".h5ad", "", data)))
    adata = read_h5ad(file.path(analysis_dir, data))

    ##
    cluster_meta_species = cluster_meta[cluster_meta$Group %in% adata$obs_names,]

    ## Arrange adata by cluster meta
    adata = adata[cluster_meta_species$Group]

    ## -------------------------------------------------
    ## Apply data normalization for plotting
    ## -------------------------------------------------

    ## min-max normalize along the columns
    adata$layers["mean_min-max"] = matrix(0, nrow = nrow(adata$layers["mean"]), ncol = ncol(adata$layers["mean"])) ## Need to initialize the matrix first
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
        tau = adata$var[[tolower(species)]],
        max_cell_type = factor(max_ct, levels = cluster_meta_species$Group),
        max_expr = max_val
    )

    # ## Arrange genes by their max cell type and keep the top 10 based on tau
    # gene_table_sorted = gene_table %>%
    #     arrange(max_cell_type, desc(tau)) %>%
    #     # group_by(max_cell_type) %>%
    #     # slice_head(n = 10) %>%
    #     # ungroup() %>%
    #     as.data.frame()

    # ## Add to adata object
    # adata$var["max_cell_type"] = gene_table_sorted$max_cell_type[match(adata$var_names, gene_table_sorted$gene)]
    # adata$var["tau"] = gene_table_sorted$tau[match(adata$var_names, gene_table_sorted$gene)]
    # adata$var["max_expr"] = gene_table_sorted$max_expr[match(adata$var_names, gene_table_sorted$gene)]
    # adata$var["gene_order"] = match(gene_table_sorted$gene, adata$var_names)

    ## Arrange genes by their max cell type and keep the top 10 based on tau
    gene_table_sorted_filtered = gene_table %>%
        arrange(max_cell_type, desc(tau)) %>%
        group_by(max_cell_type) %>%
        slice_head(n = 25) %>%
        ungroup() %>%
        as.data.frame()

    ##
    adata = adata[,gene_table_sorted_filtered$gene]

    ## Add to adata object
    adata$var["max_cell_type"] = gene_table_sorted_filtered$max_cell_type[match(adata$var_names, gene_table_sorted_filtered$gene)]
    adata$var["tau"] = gene_table_sorted_filtered$tau[match(adata$var_names, gene_table_sorted_filtered$gene)]
    adata$var["max_expr"] = gene_table_sorted_filtered$max_expr[match(adata$var_names, gene_table_sorted_filtered$gene)]
    adata$var["gene_order"] = match(gene_table_sorted_filtered$gene, adata$var_names)

    ##
    adata_set[[species]] = adata
}

## Gather union of cell types across species
cell_types = Reduce(union, lapply(adata_set, function(x) { x$obs_names } ))
cell_types = cell_types[match(cluster_meta$Group, cell_types)]  ## Reorder to match cluster meta
cell_types = factor(cell_types, levels = cluster_meta$Group)

## For any species that is missing a cell_type, add a row of zeros
for(species in names(adata_set)) {
    adata <- adata_set[[species]]
    min_max = adata$layers["mean_min-max"]
    adata = AnnData(X = min_max, obs = adata$obs, var = adata$var)
    adata$layers["mean_min-max"] = min_max
    ## Find missing cell types
    existing <- as.character(adata$obs_names)
    missing <- setdiff(cell_types, existing)
    ## If missing, add a row of zeros
    if(length(missing) > 0) {
        message("Adding missing cell types for: ", species, " -> ", paste(missing, collapse = ", "))

        n_vars <- adata$n_vars
        new_obs <- data.frame(ct = missing); rownames(new_obs) = missing
        new_X <- matrix(0, nrow=length(missing), ncol=n_vars)

        ## Create new AnnData object with missing cell types
        adata_missing <- AnnData(X = new_X, obs = new_obs, var = adata$var)

        ## concatenate along observations (rows)
        adata <- concat(list(adata, adata_missing), join = "outer", axis = 0L, merge = "first")
    }
    ## Reorder to match cell_types
    adata <- adata[match(cell_types, adata$obs_names), ]
    ## Reorder genes to match gene_table_sorted
    adata$layers["mean_min-max"] = adata$X
    adata_set[[species]] = adata
}

## Create cell type blocks and annotation bar
class_meta = cluster_meta %>%
    distinct(Class, .keep_all = TRUE)
class_colors = class_meta$Class_colors
names(class_colors) = class_meta$Class

## Now draw plots for each species gene set and across species
heatmap_list = NULL
for(species in c("Human", "Macaque", "Mouse")) {

    adata_plot = adata_set[[species]]
    adata_plot = adata_plot[,adata_plot$var$gene_order]
    adata_plot = adata_plot[,adata_plot$var_names[adata_plot$var$tau > 0.9]]

    ## Add Class info to obs
    adata_plot$obs$Class = cluster_meta$Class[match(adata_plot$obs_names, cluster_meta$Group)]
    adata_plot$var[["Class"]] = rep("", nrow(adata_plot$var))
    adata_plot$var[["Class"]] = cluster_meta$Class[match(adata_plot$var$max_cell_type, cluster_meta$Group)]

    ## Heatmap blocks
    adata_plot_neuron = adata_plot[adata_plot$obs$Class != "Non-Neurons",]
    celltype_sep = factor(adata_plot_neuron$obs$Class, levels = class_meta$Class)

    ha.celltype <- HeatmapAnnotation(
        celltype = factor(adata_plot_neuron$obs$Class, levels = class_meta$Class),
        col = list(celltype = class_colors),
        show_annotation_name = FALSE, 
        show_legend = FALSE,
        height = unit(2, "mm")
    )

    adata_plot_neuron$var %>% filter(Class != "Non-Neurons") %>% rownames() -> gene_set

    print(paste0("Drawing heatmap for species: ", species, " with ", ncol(adata_plot), " genes."))

    ## ComplexHeatmap
    heatmap = Heatmap(t(adata_plot_neuron[,adata_plot_neuron$var_names %in% gene_set]$layers["mean_min-max"]),
                    cluster_columns = FALSE, 
                    show_column_dend = FALSE,
                    cluster_rows = FALSE,
                    col=colorRampPalette(c("#FFFFFF",color_meta[color_meta$label == species,"color_hex_triplet"]))(100),
                    # top_annotation = ha.celltype, 
                    column_split = celltype_sep,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_title = NULL,
                    column_title = NULL,
                    border=TRUE,
                    na_col = "#F5F5F5",
                    show_heatmap_legend = FALSE,
                    show_column_names = if(species == "Mouse") TRUE else FALSE,
                    show_row_names=FALSE)

    ## Combine horizontally
    if (is.null(heatmap_list)) {
        heatmap_list <- heatmap
    } else {
        heatmap_list <- heatmap_list %v% heatmap
    }
}

png(paste0(file.path(figure_dir, "marker_plot_group_species_specific_neuron.png")), width=10, height=5, units="in", res=600)
draw(heatmap_list)
dev.off()

