library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(tidyr)
library(forcats)
library(anndata)

##
marker_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs"
anno_term = "group"

## Load gene sets
tf_database = read.table("/home/nelson.johansen/scripts/EvoGen/Projects/M1Evo/data/metadata/TFs_lambert_pmid29425488_1.01.txt", sep="\t", header=TRUE)
tf_database = tf_database$hgnc_symbol

##
tauResults = read.csv(file.path(marker_dir, "tau_scores", paste0(anno_term, "_tau_scores_mean_donor_meta_analysis_all.csv")))
tauResults$gene = tauResults$X; tauResults$X = NULL
expressologResults = read.csv(file.path(analysis_dir, anno_term, paste0("expressologs_", anno_term, ".csv")))

## Multi-species expressologs
expressologs = expressologResults %>% 
                group_by(gene) %>% 
                mutate(median_cor = median(expr_cor, na.rm = TRUE)) %>%
                ungroup() %>%
                arrange(desc(median_cor)) %>%
                mutate(rank = dense_rank(desc(-median_cor))) %>%
                mutate(auroc = rank / n_distinct(gene)) %>% as.data.frame()

write.csv(expressologs, file.path(analysis_dir, anno_term, paste0(anno_term, "_expressologs_with_metrics.csv")), row.names = FALSE)

## Load metadata

## AnnoTable
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet="consensus_anno_pre-print")
cluster_meta = cluster_meta %>%
    distinct(Group, .keep_all = TRUE)
cluster_meta = cluster_meta[cluster_meta$display_order_group,]
cluster_meta$groupby = cluster_meta$Group

group_hierarchy_diagram_order = cluster_meta$Group

species_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/HMBA_BG_consensus_annotation.xlsx", sheet = "value_set_colors")
species_meta$species = tolower(species_meta$label)
species_meta$species_color = species_meta$color_hex_triplet

species_meta = species_meta %>% filter(label %in% c("Human", "Macaque", "Marmoset"))

##
species_order = c("human", "macaque", "marmoset")

## -------------------------------------------
## Plot 1: Distribution of median gene correlations across groupbyes between species.
plot.me <- expressologs %>% 
  distinct(gene, .keep_all = TRUE)

pdf(paste0(analysis_dir, "/figures/expressologs_dist.pdf"), width = 8, height=8) #, units = "in", res = 300)
hist(plot.me$median_cor, breaks=25)
dev.off()

## -------------------------------------------
## Plot2:  Heatmap of mean correlation of genes (by groupby) between species.
plot.me <- expressologResults %>% 
  # filter(sp1_expr_stat > 0.2 | sp2_expr_stat > 0.2) %>% 
  group_by(species1, species2) %>% 
  summarize(mean_cor = mean(expr_cor, na.rm=T)) %>% 
  pivot_wider(names_from = species2, values_from = mean_cor) %>% 
  ungroup() %>% 
  select(-species1) %>% 
  as.matrix()

## Figure out matched species order for rows
species.idx = apply(plot.me, 2, function(x) which(is.na(x)))
rownames(plot.me) <- names(sort(species.idx))

## ComplexHeatmap without clustering
pdf(paste0(analysis_dir, "/figures/expressologs_heatmap.pdf"), width = 9, height=8) #, units = "in", res = 300)
Heatmap(plot.me[species_order, species_order], 
        col=magma(100),
        cluster_rows = FALSE, 
        cluster_columns = FALSE)
dev.off()

## -------------------------------------------
## Plot3:  Mean cor by gene/species
plot.me <- expressologResults %>% 
  filter(sp_order == 1) %>%
  # filter(sp1_expr_stat > 0.2 | sp2_expr_stat > 0.2) %>% 
  group_by(species1, gene) %>% 
  summarize(mean_cor = mean(expr_cor, na.rm=T)) %>% 
  ggplot(aes(x = factor(species1, rev(species_order)), y = mean_cor, colour=species1)) +
  geom_boxplot() +
  # scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

pdf(paste0(analysis_dir, "/figures/expressologs_mean_corr.pdf"), width = 8, height=8)
plot(plot.me)
dev.off()

## -------------------------------------------
## Plot4: Most divergent genes using human as reference
gene_order = expressologResults %>% 
  filter(species1 == "human") %>% 
  group_by(gene) %>% 
  summarize(mean_cor = mean(expr_cor, na.rm=TRUE))
gene_order = gene_order[order(-gene_order$mean_cor),]$gene

plot.me <- expressologResults %>% 
  filter(species1 == "human") %>% 
  group_by(gene) %>%
  mutate(median_cor = median(expr_cor, na.rm=TRUE)) %>%
  filter(median_cor > 0) %>%
  ggplot(aes(y = as.numeric(factor(gene, gene_order)), x = expr_cor, color = species2)) +
  scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
  #geom_point(size=0.1) +
  geom_smooth(span = 0.2, size=0.5) +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

pdf(paste0(analysis_dir, "/figures/expressologs_divergent_human.pdf"), width = 4, height=8)
plot(plot.me)
dev.off()

plot.me <- expressologResults %>% 
  filter(species1 == "human") %>% 
  group_by(gene) %>%
  mutate(median_cor = median(expr_cor, na.rm=TRUE)) %>%
  filter(median_cor > 0) %>%
  ggplot(aes(y = as.numeric(factor(gene, gene_order)), x = expr_cor, color = species2)) +
  # scale_color_manual(values = species_meta$color, breaks=species_meta$species) +
  geom_point(size=0.05) +
  geom_smooth(span = 0.2, size=0.5) +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

pdf(paste0(analysis_dir, "/figures/expressologs_divergent_human_with_points.pdf"), width = 4, height=8)
plot(plot.me)
dev.off()

## -------------------------------------------
## Plot5: Density plot of most divergent genes by species using human as reference (similar to plot4)
plot.me <- expressologResults %>% 
  filter(species1 == "human" & sp1_expr_stat > 0) %>% 
  ggplot(aes(x = expr_cor, color = species2)) +
  geom_density() +
  # scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
  theme_bw() +
  theme(legend.position = "none")
  

pdf(paste0(analysis_dir, "/figures/expressologs_divergent_human_species.pdf"), width = 4, height=8)
plot(plot.me)
dev.off()

## -------------------------------------------
## File1: Human divergent expressolog results
# human_div_genes <- expressologResults %>%
#   filter(species1 == "human") %>%
#   group_by(gene) %>% 
#   summarize(mean_cor = mean(expr_cor)) %>% 
#   arrange(mean_cor)

# write_csv(human_div_genes, file = "Projects/M1evo/M1evo-RNA/paper_analysis/expressologs/expressologs_human_div_gene.csv")

## -------------------------------------------
## Data1: Setup data for expressolog line plots

expr_summary = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/mean_expr/group"
expr_files = list.files(expr_summary, pattern="*_group_mean_expr_orthologs.h5ad")

## Load expression data from hdf
data = list()
for(file in expr_files){
  adata = read_h5ad(file.path(expr_summary, file))
  species_name = gsub("_group_mean_expr_orthologs\\.h5ad$", "", file)
  data[[species_name]] = t(as.data.frame(adata$layers[["mean"]]))
}

## -------------------------------------------
## Plot6: Expressolog line plots for a gene of interest

#' Expressolog line plots for a gene of interest
#'
#' @param gene.interest Gene of interest.
#' @param data List of expression data for each species.
#' @param species_meta Species metadata.
#' @param cluster_meta Cluster metadata.
#' @param expressologs Expressologs data.
#'
#' @return Line plots for gene of interest.
expresso_lines = function(gene.interest, data, species_meta, cluster_meta, expressologs.data, plot.width=10, plot.height=3){

  ##
  df.species = NULL
  for(species in unique(expressologs.data$species1)){
      df.human = as.data.frame(data[[species]])
      df.human$genes = gsub("GENE-", "", toupper(rownames(df.human)))
      if(gene.interest %in% df.human$genes){
        df.human <- df.human %>%
          filter(genes == !!gene.interest) %>%
          pivot_longer(cols = -genes, names_to = "groupby", values_to = "expr")
        df.human$species=species
      }else{
        df.human <- data.frame(genes = gene.interest, groupby = colnames(df.human), expr = NA, species = species)
      }
      df.species = rbind(df.species, df.human)
  }

  ## Compute AUROC on the fly
  gene_AUROC = expressologs.data %>% filter(gene == gene.interest) %>% pull(rank) %>% first() / length(unique(expressologs.data$gene))

  ## Plot lines coloring by species
  plot.me <- df.species %>% 
      inner_join(cluster_meta, by = "groupby") %>%
      inner_join(species_meta, by = "species") %>%
      mutate(groupby = fct_reorder(groupby, display_order_group)) %>%
      # slice_min(mean_cor, n = 12*22) %>%
      filter(genes == !!gene.interest) %>% 
      mutate(expr = replace_na(expr, 0)) %>%
      # mutate(gene = fct_reorder(gene, mean_cor)) %>% 
      # ggplot(aes(x = groupby, y = expr, colour = phylo1_sci, group=species)) +
      ggplot(aes(x = groupby, y = expr, colour = species, group = species)) +
      scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
      geom_point() +
      geom_line(aes(linetype=species)) +
      ggtitle(paste0(gene.interest, ": ", round(gene_AUROC, 4), "(Conservation (auroc))")) +
      xlab("Group") +
      ylab("Expression") +
      theme_bw() +
      theme(axis.text = element_text(color = "#000000"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              panel.background = element_rect(fill = "#FFFFFF"))

  pdf(paste0(analysis_dir, "/figures/expressolog_lines/expressologs_lines_species_", gene.interest, ".pdf"), width=plot.width, height=plot.height)
  plot(plot.me)
  dev.off()

  # ## Plot lines coloring by phylo1_sci
  # plot.me <- df.species %>% 
  #     inner_join(cluster_meta, by = "groupby") %>%
  #     inner_join(species_meta, by = "species") %>%
  #     mutate(groupby = fct_reorder(groupby, groupby_id)) %>%
  #     # slice_min(mean_cor, n = 12*22) %>%
  #     filter(genes == !!gene.interest) %>% 
  #     mutate(expr = replace_na(expr, 0)) %>%
  #     # mutate(gene = fct_reorder(gene, mean_cor)) %>% 
  #     ggplot(aes(x = groupby, y = expr, colour = phylo1_sci, group=species)) +
  #     scale_color_manual(values = species_meta$phylo1_color, breaks=species_meta$phylo1_sci) +
  #     geom_point() +
  #     geom_line(aes(linetype=phylo1_sci)) +
  #     ggtitle(paste0(gene.interest, ": ", round(gene_AUROC, 4), "(AUROC)")) +
  #     xlab("Subclass") +
  #     ylab("Expression") +
  #     theme_bw() +
  #     theme(axis.text = element_text(color = "#000000"),
  #             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #             panel.border = element_blank(), 
  #             panel.grid.major = element_blank(),
  #             panel.grid.minor = element_blank(), 
  #             axis.line = element_line(colour = "black"),
  #             panel.background = element_rect(fill = "#FFFFFF"))

  # pdf(paste0("Projects/M1evo/M1evo-RNA/paper_analysis/expressologs/figures/expressologs/expressologs_lines_phylo1_", gene.interest, ".pdf"), width=plot.width, height=plot.height)
  # plot(plot.me)
  # dev.off()
}


for(gene in unique(expressologs$gene)){
  print(gene)
  ## Check if file exists first
  # if(file.exists(paste0(analysis_dir, "/figures/expressolog_lines/expressologs_lines_species_", gene, ".pdf"))){
  #   next
  # }
  expresso_lines(gene, data, species_meta, cluster_meta, expressologs, plot.height=5)
}

## ---------------------------------------
## Plot7: Dot plot of expression for a set of genes in a given groupby

gene.set = tf_database
celltype.interest = "L5 ET"

## Pull data of interest and account for missing genes in certian species
# data.species = NULL
# for(species in species_meta$species_alt2){
#     medianExpr = data[[species]]
#     rownames(medianExpr) = gsub("GENE-", "", toupper(rownames(medianExpr)))
#     for(gene in gene.set){
#       if(!gene %in% rownames(medianExpr)){
#         ## Create 1 x cols matrix for new gene with NAs
#         gene.add=matrix(rep(NA, ncol(medianExpr)), nrow=1, ncol=ncol(medianExpr)); rownames(gene.add) = gene
#         medianExpr = rbind(medianExpr, gene.add)
#       }
#     }
#     data.species[[species]] = medianExpr[gene.set,celltype.interest]
# }
# data.plot = Reduce(cbind, data.species)
# colnames(data.plot) = names(data.species)

# ##
# data.plot = setNames(melt(data.plot), c('genes', 'species', 'expr'))

# ##
# plot.me = data.plot %>% 
#   filter(grepl("POU", data.plot$genes)) %>%
#   ggplot(aes(x = factor(species, species_meta$species), y = genes, color = expr, size = expr)) + 
#   geom_point() + 
#   scale_color_viridis_c(name = 'Expression') + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ylab('') +
#   theme(axis.ticks = element_blank()) 

# pdf("Projects/M1evo/M1evo-RNA/paper_analysis/expressologs/figures/TF_expression_dots.pdf", width = 10, height=10)
# plot(plot.me)
# dev.off()

# plot.me <- expressologs %>% 
#   filter(grepl("^POU", gene)) %>%
#   distinct(gene, .keep_all = TRUE)

# pdf("Projects/M1evo/M1evo-RNA/paper_analysis/expressologs/figures/expressologs_dist_TF.pdf", width = 8, height=8) #, units = "in", res = 300)
# Heatmap(as.matrix(plot.me$auroc[rev(order(plot.me$gene))]), col=magma(100), cluster_rows = FALSE, cluster_columns = FALSE)
# dev.off()



# matrix_data <- data.plot %>%
#   pivot_wider(names_from = species, values_from = expr) %>%
#   column_to_rownames(var = "genes")

# na.counts = colSums(apply(matrix_data,1,is.na))
# na.counts = names(na.counts[na.counts < 10])
# matrix_data = matrix_data[na.counts,]

all_columns = unique(unlist(lapply(data, colnames)))
all_rows = unique(unlist(lapply(data, rownames)))

matrix_data = matrix(0, nrow=length(all_rows), ncol=length(all_columns), dimnames=list(all_rows, all_columns))
for(gene in all_rows){
  for(celltype in all_columns){
    medianExpr = c()
    for(species in names(data)){
      if(celltype %in% colnames(data[[species]]) && gene %in% rownames(data[[species]])){
        medianExpr = c(medianExpr, data[[species]][gene, celltype])
      }
      if(!is.na(medianExpr) && length(medianExpr) > 0){
        ## Compute median expression across species
        medianExpr = median(medianExpr, na.rm = TRUE)
      }else{
        medianExpr = 0
      }
      matrix_data[gene, celltype] = median(medianExpr)
    }
  }
}

##
auroc_tau = merge(tauResults, expressologs, by = "gene", all.x = TRUE) 
sub_expressologs = auroc_tau %>% filter(auroc > 0.75) %>% filter(xspecies_min > 0.75) %>% pull(gene)
sub_expressologs = unique(sub_expressologs)


##
sub_matrix = matrix_data[sub_expressologs,unique(cluster_meta$Group)]

##
cluster_meta_unique <- cluster_meta %>% distinct(Group, .keep_all = TRUE)

# hclust_data = matrix_data
# hclust_data[is.na(hclust_data)] = 0
# ordering = hclust(dist(hclust_data))$order

png(paste0(analysis_dir, "/figures/TF_expression_heatmap.png"), units="in", res=600, width = 10, height=10)
Heatmap(sub_matrix,
        col=magma(100),
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 4),
        column_split = factor(cluster_meta_unique$Class, levels = c("Motor Neurons", "Cholinergic", "Glutamatergic", "GABAergic", "GABA-Gluta", "Non-Neurons")),
        na_col = "gray")
dev.off()


## matrix_data[c(names(sort(-apply(matrix_data,1,var)))[1:100], "POU3F1"),]