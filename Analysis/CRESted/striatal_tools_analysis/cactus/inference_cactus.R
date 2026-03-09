library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(anndata)
library(openxlsx)
library(viridisLite)
library(scales)

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

striatal_table = read.xlsx(file.path(analysis_dir, "crested", "enhancer_tables", "striatal_tools_supp_table1.xlsx"))

## -------------------------
## Load data
## -------------------------
adata = read_h5ad(file.path(analysis_dir, "crested", "data", "striatal_tools_enhancer_maf_predictions.h5ad"))
# cactus_adata = read_h5ad(file.path(analysis_dir, "cactus", "analysis", "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"), backed="r")

# ## Species order 
# species_order = as.character(cactus_adata$var$primates_shorthand)
# species_order[is.na(species_order)] = as.character(cactus_adata$var$order[is.na(species_order)])
# species_order[is.na(species_order)] = "unknown"

species_meta = read.xlsx(file.path(analysis_dir, "cactus", "CACTUS_species_info.xlsx"))
species_clade = read.xlsx(file.path(analysis_dir, "cactus", "CACTUS_species_info.xlsx"), sheet="clade_order")

## Define levels for species order
species_order = as.character(species_meta$primates_shorthand)
species_order[is.na(species_order)] = as.character(species_meta$order[is.na(species_order)])

## -------------------------
## Define plotting orders
## -------------------------
species_order = factor(species_order, levels=species_clade$clade)

adata$var["clade"] = species_order

## -------------------------
## Plot Heatmap
## -------------------------

# enh_order = c("AiE0612h", "AiE0621h", "AiE0618h","AiE0951h", "AiE0445h", "AiE0351h", "AiE0617h",
#                 "AiE0779m", "AiE0603m", "AiE1190m", "AiE0868m", "AiE0600m", "AiE0785m",
#                 "AiE0873m", "AiE1404m", "AiE0869m",
#                 "AiE0441h", "AiE0452h", "AiE0447h",
#                 "AiE0780m", "AiE0782m", "AiE0784m",
#                 "AiE0450h", "AiE0616h", 
#                 "AiE0769m", "AiE1191m",
#                 "AiE0444h", "AiE1026h", "AiE0367h",
#                 "AiE1421m", "AiE1426m",
#                 "AiE0682h", 
#                 "AiE1419m", "AiE0743m")

## Extract species from end of enhancer name, e.g. AiE0001h -> h (human) and AiE0001m -> m (mouse)
enhancer_species = sapply(adata$obs_names, function(x) {
    strsplit(x, split="")[[1]][nchar(x)]
})
adata$obs["species"] = enhancer_species
adata$obs["species"] = gsub("m", "mouse", adata$obs$species)
adata$obs["species"] = gsub("h", "human", adata$obs$species)

## 
enhancer_celltype = sapply(adata$obs_names, function(x) {
    striatal_table[which(gsub(" ", "", striatal_table[["enhancer.ID"]]) == x), "Cell.type"]
})
adata$obs["enhancer_celltype"] = enhancer_celltype
adata$obs["enhancer_celltype"] = factor(adata$obs[["enhancer_celltype"]], levels=c("pan-MSN", "D1 MSN", "D2 MSN", "Cholinergic", "Sst-Chodl", "Pvalb-Pthlh"))

## order adata by species
for(species in c("mouse", "human")) {
    plot_adata = adata[which(adata$obs$species == species),]$copy() #[match(enh_order, adata$obs_names),]$copy()

    ## Annotation
    ha.species <- rowAnnotation(
        species = plot_adata$obs$species,
        enhancer = plot_adata$obs$celltype,
        show_annotation_name = FALSE, 
        show_legend = FALSE
    )

    ## Normalize data matrix
    norm_matrix = apply(plot_adata$X, 1, function(x) {
        (x - min(x)) / ((max(x) - min(x)) + 1e-5)
    })
    plot_adata$layers["min-max"] = t(norm_matrix)

    ## Replace 0 with NA for better visualization
    plot_adata$layers["min-max"][plot_adata$layers["min-max"] == 0] = NA

    ## Colors
    # breaks <- seq(0, 1, length.out = 100)
    ## Blue to Red color map
    # col_fun <- colorRamp2(c(0, 1), c("blue", "red"))

    breaks <- seq(0, 1, length.out = 100)
    col_fun <- colorRamp2(breaks, magma(length(breaks)))

    ## ComplexHeatmap Neurons
    heatmap = Heatmap(plot_adata$layers["min-max"],
                    cluster_columns = FALSE, 
                    show_column_dend = FALSE,
                    cluster_rows = FALSE,
                    show_row_dend = FALSE,
                    col=col_fun,
                    right_annotation = ha.species,
                    # bottom_annotation = ha.tf_hit,
                    # top_annotation = ha.celltype, 
                    row_split=plot_adata$obs["enhancer_celltype"],
                    column_split=species_order,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    na_col = "grey80",
                    row_title = NULL,
                    column_title = NULL,
                    border=TRUE,
                    show_heatmap_legend = FALSE,
                    show_column_names = TRUE,
                    show_row_names= TRUE)

    pdf(file.path(analysis_dir, "crested", "figures", paste0("striatal_tools_", species, "_enhancers_maf_predictions.pdf")), width=20, height=8)
    draw(heatmap)
    dev.off()
}

## -------------------------
## Main figure 
## -------------------------
plot_adata = adata[which(adata$obs_names == "AiE0779m"),]$copy()

## Normalize values between 0 and 1
values = as.numeric(plot_adata$X[1,])
min_val <- min(values, na.rm = TRUE)
max_val <- max(values, na.rm = TRUE)
values_norm <- (values - min_val) / (max_val - min_val)

## Aggregate (mean) per clade
df <- data.frame(
  clade = adata$var["clade"],
  value = values_norm
)

## Calculate mean per clade
clade_means <- df %>%
                    group_by(clade) %>%
                    summarise(mean_value = mean(value, na.rm = TRUE)) %>%
                    as.data.frame()

## Plot
df$shape <- ifelse(df$value > 0.5, 22, 21)

# You can also map size to value
max_point_size <- 8
df$point_size <- rescale(df$value, to = c(2, max_point_size))

# Define color gradient
col_low <- "#fde0dd"
col_high <- "#de2d26"

clades <- unique(df$clade)
rect_df <- data.frame(
  ymin = seq_along(clades) - 0.5, 
  ymax = seq_along(clades) + 0.5,
  xmin = -Inf,                    
  xmax = Inf
)

ggplot(df, aes(x = value, y = factor(clade, levels = rev(levels(factor(clade)))))) +
  geom_rect(
    data = rect_df,
    aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
    inherit.aes = FALSE,
    fill = "white",
    alpha = 0.05  # transparent bands
  ) +
  geom_hline(
    yintercept = seq_along(clades) + 0.5,
    color = "gray60",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(size = point_size, fill = value),
    shape = 21,
    color = "black",
    stroke = 0.4,
    position = position_jitter(width = 0.25, height = 0.25)
  ) +
  scale_shape_identity() +
  scale_fill_viridis_c(option = "plasma") +
  scale_size_identity() +
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  guides(size = "none", shape = "none") +
  labs(fill = "Value")

ggsave(file.path(analysis_dir, "crested", "figures", "striatal_tools_AiE0779m_enhancer_maf_predictions_summary.pdf"), width=6, height=12, units="in", dpi=900)

## Identify max per clade
df$species <- rownames(df)

# Get max row per clade
clade_max <- df %>%
  group_by(clade) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()