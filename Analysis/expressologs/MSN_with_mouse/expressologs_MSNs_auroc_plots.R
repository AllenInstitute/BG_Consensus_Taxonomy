library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(ggplot2)
library(ggrepel)

##
marker_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/marker_gene_analysis/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs/"
anno_term = "group"

## Load gene sets
tf_database = read.table("/home/nelson.johansen/scripts/EvoGen/Projects/M1Evo/data/metadata/TFs_lambert_pmid29425488_1.01.txt", sep="\t", header=TRUE)
tf_database = tf_database$hgnc_symbol

##
tauResults = read.csv(file.path(marker_dir, "tau_scores", "MSN_vignette", paste0(anno_term, "_tau_scores_mean_donor_meta_analysis_all.csv")))
tauResults$gene = tauResults$X; tauResults$X = NULL

expressologResults = read.csv(file.path(analysis_dir, "MSN_vignette", anno_term, "expressologs_MSN_group.csv"))

## Multi-species expressologs
expressologs = expressologResults %>% 
                group_by(gene) %>% 
                mutate(median_cor = median(expr_cor, na.rm = TRUE)) %>%
                ungroup() %>%
                arrange(desc(median_cor)) %>%
                mutate(rank = dense_rank(desc(-median_cor))) %>%
                mutate(auroc = rank / n_distinct(gene)) %>% as.data.frame()

## Merge tau and expressologs
auroc_tau = merge(tauResults, expressologs, by = "gene", all.x = TRUE) 

##
auroc_tau = auroc_tau %>% filter(!is.na(auroc) & !is.na(xspecies_min))

## Save
write.csv(auroc_tau, file.path(analysis_dir, "MSN_vignette", anno_term, paste0(anno_term, "_expressologs_MSN_auroc_tau_metrics.csv")), row.names = FALSE)

##
plot_two_pairs_without_geneexpr <- function(final_plot_df, x_val, y_val, x_axis_top, y_axis_top, TF_database=NULL) {

  # Get top 10 genes by y_val
  top_genes <- final_plot_df %>%
              # Arrange by both x_val and y_val in descending order
              arrange(desc(.data[[x_val]]), desc(.data[[y_val]])) %>%
              # Take the top 10 rows
              slice(1:10) %>%
              # Remove duplicates based on gene, keep all other columns intact
              distinct(gene, .keep_all = TRUE)

  plot = ggplot(final_plot_df, aes(x = .data[[x_val]], y = .data[[y_val]])) +
    geom_rect(xmin=0, xmax=0.75, ymin=0, ymax=0.75, fill="#ffe6e6", alpha=0.5) + # bottom-left 
    geom_rect(xmin=0.75, xmax=1, ymin=0.75, ymax=1, fill="#e6ffe6", alpha=0.5) + # top-right "#e6ffe6"
    geom_rect(xmin=0, xmax=0.75, ymin=0.75, ymax=1, fill="#e6f7ff", alpha=0.5) + # top-left 
    geom_rect(xmin=0.75, xmax=1, ymin=0, ymax=0.75, fill="#fff5e6", alpha=0.5) + # bottom-right
    ## Now points
    geom_point(alpha = 0.5, size = 0.3, color = "#4dac26") + 
    # geom_smooth(method = "loess", se = FALSE, color = "blue", linetype = "dashed") +
    # geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = Inf) +
    geom_vline(xintercept = x_axis_top, linetype = "dashed", color = "black") +   # vertical line
    geom_hline(yintercept = y_axis_top, linetype = "dashed", color = "black") +   # horizontal line
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Conservation (AUROC)",
         y = "Cell type specificity (Tau)") +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text = element_text(color = "#000000"),
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.background = element_rect(fill = "#FFFFFF"))
    #coord_fixed(ratio = 1)
  return(plot)
}

## Filter
auroc_tau = auroc_tau[!duplicated(auroc_tau$gene),]
auroc_tau = auroc_tau[auroc_tau$xspecies_min > 0,]

## TF subset
TF_auroc_tau = auroc_tau[auroc_tau$gene %in% tf_database,]

## Save plots
plot.me = plot_two_pairs_without_geneexpr(auroc_tau, "auroc", "xspecies_min", 0.75, 0.75)
plot.me = plot.me + geom_point(data = TF_auroc_tau, aes(x = auroc, y = xspecies_min), color = "black", fill="#d01c8b", shape=21, stroke=0.05, size = 1.5)

##
ggsave(paste0(analysis_dir, "/figures/expressologs_MSN_auroc_", "xspecies_min", "_scatterplot_final.png"), plot = plot.me, width = 4, height = 6, dpi = 900)

## Figure out the number of genes and TFs in each quadrant
n_genes_bottom_left = nrow(auroc_tau[auroc_tau$auroc <= 0.75 & auroc_tau$xspecies_min <= 0.75,])
n_genes_top_right = nrow(auroc_tau[auroc_tau$auroc > 0.75 & auroc_tau$xspecies_min > 0.75,])
n_genes_top_left = nrow(auroc_tau[auroc_tau$auroc <= 0.75 & auroc_tau$xspecies_min > 0.75,])
n_genes_bottom_right = nrow(auroc_tau[auroc_tau$auroc > 0.75 & auroc_tau$xspecies_min <= 0.75,])

## Print
print(paste0("Number of genes in bottom-left quadrant: ", n_genes_bottom_left))
print(paste0("Number of genes in top-right quadrant: ", n_genes_top_right))
print(paste0("Number of genes in top-left quadrant: ", n_genes_top_left))
print(paste0("Number of genes in bottom-right quadrant: ", n_genes_bottom_right))

## TFs
n_TFs_bottom_left = nrow(TF_auroc_tau[TF_auroc_tau$auroc <= 0.75 & TF_auroc_tau$xspecies_min <= 0.75,])
n_TFs_top_right = nrow(TF_auroc_tau[TF_auroc_tau$auroc > 0.75 & TF_auroc_tau$xspecies_min > 0.75,])
n_TFs_top_left = nrow(TF_auroc_tau[TF_auroc_tau$auroc <= 0.75 & TF_auroc_tau$xspecies_min > 0.75,])
n_TFs_bottom_right = nrow(TF_auroc_tau[TF_auroc_tau$auroc > 0.75 & TF_auroc_tau$xspecies_min <= 0.75,])

## Print
print(paste0("Number of TFs in bottom-left quadrant: ", n_TFs_bottom_left))
print(paste0("Number of TFs in top-right quadrant: ", n_TFs_top_right))
print(paste0("Number of TFs in top-left quadrant: ", n_TFs_top_left))
print(paste0("Number of TFs in bottom-right quadrant: ", n_TFs_bottom_right))