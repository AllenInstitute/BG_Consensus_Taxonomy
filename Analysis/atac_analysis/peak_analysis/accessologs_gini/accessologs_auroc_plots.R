library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(ggplot2)
library(ggrepel)

## Set directories
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Accessologs/"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/annotations"
anno_term = "Group"

## Load TF database
tf_database = read.table("/home/nelson.johansen/scripts/EvoGen/Projects/M1Evo/data/metadata/TFs_lambert_pmid29425488_1.01.txt", sep="\t", header=TRUE)
tf_database = tf_database$hgnc_symbol

## Load annotations
promoter_calls = read.csv(file.path(anno_dir, "human_peaks_promoter.csv"))
TE_calls = read.csv(file.path(anno_dir, "human_peaks_annotated_TE.csv"))
phyloP_calls = read.csv(file.path(anno_dir, "human_peaks_annotated_phyloP.csv"))

## Load ortholog calls
ortholog_calls = read.csv(file.path(anno_dir, "human_ref_liftover_HALPER_minMatch_0-5.tsv") , sep="\t")
ortholog_calls = ortholog_calls %>% filter(ortholog == "True")

## Load ATAC-seq data
spec_atac = read.csv(file.path(base_dir, "specificity", "gini_scores_combined.csv"))
                
## Replace region names with ortholog (human) names
spec_atac$human_ID = NA
for(species in unique(spec_atac$species)) {
  print(species)
  spec_atac$human_ID[spec_atac$species == species] = ortholog_calls$human_ID[match(spec_atac$region[spec_atac$species == species], ortholog_calls[[paste0(species, "_ID")]])]
}
spec_atac = spec_atac %>% filter(human_ID %in% ortholog_calls$human_ID) 

## For each unique human_ID, compute the min gini score
spec_atac = spec_atac %>%
              group_by(human_ID) %>%
              mutate(gini_scores = min(gini_scores, na.rm = TRUE)) %>%
              ungroup() %>%
              arrange(gini_scores) %>%
              distinct(human_ID, .keep_all = TRUE) %>%
              as.data.frame()

spec_atac$species_region = spec_atac$region
spec_atac$region = spec_atac$human_ID

##
accessologResults = read.csv(file.path(analysis_dir, "accessologs_mammalian_orthologs.csv"))

## Multi-species accessologs
accessologs = accessologResults %>% 
                group_by(region) %>% 
                mutate(median_cor_R = median(expr_cor, na.rm = TRUE)) %>%
                ungroup() %>%
                arrange(desc(median_cor_R)) %>%
                mutate(rank_R = dense_rank(desc(-median_cor_R))) %>%
                mutate(auroc_R = rank_R / n_distinct(region)) %>% as.data.frame()

## Merge tau and accessologs
auroc_tau = merge(spec_atac, accessologs, by = "region") 
auroc_tau = auroc_tau %>% distinct(region, .keep_all = TRUE) %>% as.data.frame()

## Merge promotor calls
auroc_tau = merge(auroc_tau, promoter_calls[,c("Name", "promoter")], by.x = "region", by.y = "Name", all.x = TRUE)
auroc_tau = merge(auroc_tau, TE_calls[,c("Name", "TE", "repName", "repClass", "repFamily")], by.x = "region", by.y = "Name", all.x = TRUE)
auroc_tau = merge(auroc_tau, phyloP_calls[,c("Name", "phyloP_mean")], by.x = "region", by.y = "Name", all.x = TRUE)

##
plot_two_pairs_without_expr <- function(final_plot_df, x_val, y_val, x_axis_top, y_axis_top, TF_database=NULL) {

  plot = ggplot(final_plot_df, aes(x = .data[[x_val]], y = .data[[y_val]])) +
    geom_point(alpha = 0.5, size = 0.5, color = "#4dac26") + # all points
    geom_vline(xintercept = x_axis_top, linetype = "dashed", color = "black") +   # vertical line
    geom_hline(yintercept = y_axis_top, linetype = "dashed", color = "black") +   # horizontal line
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Conservation (AUROC)",
         y = "Cell type specificity (Gini score)") +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text = element_text(color = "#000000"),
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.background = element_rect(fill = "#FFFFFF")) +
    coord_fixed(ratio = 1)
  return(plot)
}

##
plot_auroc_tau = auroc_tau
plot_auroc_tau <- plot_auroc_tau[!duplicated(plot_auroc_tau[c("region", "auroc")]), ]

## Promoter subset
promoter_auroc_tau = plot_auroc_tau[plot_auroc_tau$promoter == "True",]

## TE subset
TE_auroc_tau = plot_auroc_tau[plot_auroc_tau$TE == "True",]

## Save plots
plot.me = plot_two_pairs_without_expr(plot_auroc_tau, "auroc", "gini_scores", 0.75, 0.75)
plot.me = plot.me + 
            geom_point(data = promoter_auroc_tau, aes(x = auroc, y = gini_scores), color = "black", fill="#d01c8b", shape=21, stroke=0.05, size = 1.5)
plot.me = plot.me + 
            geom_point(data = TE_auroc_tau, aes(x = auroc, y = gini_scores), color = "black", fill="#f1b6da", shape=24, stroke=0.05, size = 1.5) +
            guides(fill = guide_legend(override.aes = list(size=5))) +
            theme(legend.position = "right") +
            labs(title = "Ortholog cCRE ATAC-seq: Accessologs AUROC vs Gini scores",
                 subtitle = paste0("n = ", nrow(plot_auroc_tau), " cCREs; Promoter = ", nrow(promoter_auroc_tau), " cCREs; TE = ", nrow(TE_auroc_tau), " cCREs"),
                 caption = "Promoter = pink circles; TE = light pink triangles")  

##
ggsave(paste0(analysis_dir, "/figures/accessologs_auroc_", "gini_scores_macaque", "_scatterplot.png"), plot = plot.me, width = 16, height = 12, dpi = 600)

## Save data
write.csv(auroc_tau, paste0(analysis_dir, "/accessologs_with_human_gini_annotations.csv"), row.names = FALSE)
