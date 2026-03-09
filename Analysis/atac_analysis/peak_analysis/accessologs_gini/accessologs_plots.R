library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(ggplot2)
library(ggrepel)
library(openxlsx)

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Accessologs/"
anno_term = "Group"

##
spec_atac = read.csv(file.path(base_dir, "gini_scores_combined.csv"))
spec_atac$X = NULL
spec_atac = spec_atac %>% 
                filter(species == "human")

##
accessologResults = read.csv(file.path(analysis_dir, "accessologs_mammalian_orthologs.csv"))

## Multi-species accessologs
accessologs = accessologResults %>% 
                group_by(gene) %>% 
                mutate(median_cor_R = median(expr_cor, na.rm = TRUE)) %>%
                ungroup() %>%
                arrange(desc(median_cor_R)) %>%
                mutate(rank_R = dense_rank(desc(-median_cor_R))) %>%
                mutate(auroc_R = rank_R / n_distinct(gene)) %>% as.data.frame()

## Load metadata
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta$groupby = cluster_meta$Group
species_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx", sheet = "species_metadata")
species_meta$species = tolower(species_meta$species)

##
species_order = c("human", "macaque", "mouse")

## -------------------------------------------
## Plot 1: Distribution of median gene correlations across groupbyes between species.
plot.me <- accessologs %>% 
  distinct(gene, .keep_all = TRUE)

pdf(paste0(analysis_dir, "figures/accessologs_dist.pdf"), width = 8, height=8) #, units = "in", res = 300)
hist(plot.me$median_cor, breaks=25)
dev.off()

## -------------------------------------------
## Plot2:  Heatmap of mean correlation of genes (by groupby) between species.
plot.me <- accessologResults %>% 
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
pdf(paste0(analysis_dir, "/figures/accessologs_heatmap.pdf"), width = 9, height=8) #, units = "in", res = 300)
Heatmap(plot.me[species_order, species_order], 
        col=magma(100),
        cluster_rows = FALSE, 
        cluster_columns = FALSE)
dev.off()

## -------------------------------------------
## Plot3:  Mean cor by gene/species
plot.me <- accessologResults %>% 
  filter(sp_order == 1) %>%
  # filter(sp1_expr_stat > 0.2 | sp2_expr_stat > 0.2) %>% 
  group_by(species1, gene) %>% 
  summarize(mean_cor = mean(expr_cor, na.rm=T)) %>% 
  ggplot(aes(x = factor(species1, rev(species_order)), y = mean_cor, colour=species1)) +
  geom_boxplot() +
  scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

pdf(paste0(analysis_dir, "/figures/accessologs_mean_corr.pdf"), width = 8, height=8)
plot(plot.me)
dev.off()

## -------------------------------------------
## Plot4: Most divergent genes using human as reference
gene_order = accessologResults %>% 
  filter(species1 == "human") %>% 
  group_by(gene) %>% 
  summarize(mean_cor = mean(expr_cor, na.rm=TRUE))
gene_order = gene_order[order(-gene_order$mean_cor),]$gene

plot.me <- accessologResults %>% 
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

pdf(paste0(analysis_dir, "/figures/accessologs_divergent_human.pdf"), width = 4, height=8)
plot(plot.me)
dev.off()

plot.me <- accessologResults %>% 
  filter(species1 == "human") %>% 
  group_by(gene) %>%
  mutate(median_cor = median(expr_cor, na.rm=TRUE)) %>%
  filter(median_cor > 0) %>%
  ggplot(aes(y = as.numeric(factor(gene, gene_order)), x = expr_cor, color = species2)) +
  scale_color_manual(values = species_meta$color, breaks=species_meta$species) +
  geom_point(size=0.05) +
  geom_smooth(span = 0.2, size=0.5) +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

png(paste0(analysis_dir, "/figures/accessologs_divergent_human_with_points.png"), width = 4, height=12, units="in", res=600)
plot(plot.me)
dev.off()

## -------------------------------------------
## Plot5: Density plot of most divergent genes by species using human as reference (similar to plot4)
plot.me <- accessologResults %>% 
  filter(species1 == "human" & sp1_expr_stat > 0) %>% 
  ggplot(aes(x = expr_cor, color = species2)) +
  geom_density() +
  # scale_color_manual(values = species_meta$species_color, breaks=species_meta$species) +
  theme_bw() +
  theme(legend.position = "none")
  

pdf(paste0(analysis_dir, "/figures/accessologs_divergent_human_species.pdf"), width = 4, height=8)
plot(plot.me)
dev.off()

## -------------------------------------------
## File1: Human divergent accessolog results
# human_div_genes <- accessologResults %>%
#   filter(species1 == "human") %>%
#   group_by(gene) %>% 
#   summarize(mean_cor = mean(expr_cor)) %>% 
#   arrange(mean_cor)

# write_csv(human_div_genes, file = "Projects/M1evo/M1evo-RNA/paper_analysis/accessologs/accessologs_human_div_gene.csv")

## -------------------------------------------
## Data1: Setup data for accessolog line plots

expr_summary = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/accessologs/"
expr_files = list.files(expr_summary, pattern="*_group_mean_expr.h5ad")

## Load expression data from hdf
data = list()
for(file in expr_files){
  adata = read_h5ad(file.path(expr_summary, file))
  species_name = gsub("_group_mean_expr\\.h5ad$", "", file)
  data[[species_name]] = t(as.data.frame(adata$layers[["mean"]]))
}

## -------------------------------------------
## Plot6: accessolog line plots for a gene of interest

#' accessolog line plots for a gene of interest
#'
#' @param gene.interest Gene of interest.
#' @param data List of expression data for each species.
#' @param species_meta Species metadata.
#' @param cluster_meta Cluster metadata.
#' @param accessologs accessologs data.
#'
#' @return Line plots for gene of interest.
expresso_lines = function(gene.interest, data, species_meta, cluster_meta, accessologs.data, plot.width=10, plot.height=3){

  ##
  df.species = NULL
  for(species in unique(accessologs.data$species1)){
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
  gene_AUROC = accessologs.data %>% filter(gene == gene.interest) %>% pull(rank) %>% first() / length(unique(accessologs.data$gene))

  ## Plot lines coloring by species
  plot.me <- df.species %>% 
      inner_join(cluster_meta, by = "groupby") %>%
      inner_join(species_meta, by = "species") %>%
      mutate(groupby = fct_reorder(groupby, group_hierarchy_diagram_order)) %>%
      # slice_min(mean_cor, n = 12*22) %>%
      filter(genes == !!gene.interest) %>% 
      mutate(expr = replace_na(expr, 0)) %>%
      # mutate(gene = fct_reorder(gene, mean_cor)) %>% 
      # ggplot(aes(x = groupby, y = expr, colour = phylo1_sci, group=species)) +
      ggplot(aes(x = groupby, y = expr, colour = species, group = species)) +
      scale_color_manual(values = species_meta$color, breaks=species_meta$species) +
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

  png(paste0(analysis_dir, "/figures/accessologs_TF/accessologs_lines_species_", gene.interest, ".png"), res=600, units="in", width=plot.width, height=plot.height)
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

  # pdf(paste0("Projects/M1evo/M1evo-RNA/paper_analysis/accessologs/figures/accessologs/accessologs_lines_phylo1_", gene.interest, ".pdf"), width=plot.width, height=plot.height)
  # plot(plot.me)
  # dev.off()
}


for(gene in unique(tf_database)){
  print(gene)
  ## Check if file exists first
  if(file.exists(paste0(analysis_dir, "/figures/accessologs_TF/accessologs_lines_species_", gene, ".png"))){
    next
  }
  expresso_lines(gene, data, species_meta, cluster_meta, accessologs)
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

# pdf("Projects/M1evo/M1evo-RNA/paper_analysis/accessologs/figures/TF_expression_dots.pdf", width = 10, height=10)
# plot(plot.me)
# dev.off()

# plot.me <- accessologs %>% 
#   filter(grepl("^POU", gene)) %>%
#   distinct(gene, .keep_all = TRUE)

# pdf("Projects/M1evo/M1evo-RNA/paper_analysis/accessologs/figures/accessologs_dist_TF.pdf", width = 8, height=8) #, units = "in", res = 300)
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
auroc_tau = merge(tauResults, accessologs, by = "gene", all.x = TRUE) 
sub_accessologs = auroc_tau %>% filter(auroc > 0.75) %>% filter(xspecies_min > 0.75) %>% pull(gene)
sub_accessologs = unique(sub_accessologs)


##
sub_matrix = matrix_data[sub_accessologs,unique(cluster_meta$Group)]

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