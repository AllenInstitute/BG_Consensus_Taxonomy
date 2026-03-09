library(ggalluvial)
library(openxlsx)
library(dplyr)
library(ggplot2)

analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/homology/figures/"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## -------------------------
## AnnoTable
## -------------------------
cluster_meta = read.xlsx(file.path(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet="consensus_anno_pre-print")
cluster_meta = cluster_meta %>% distinct(Group, .keep_all=TRUE)

##
homology = read.xlsx(file.path(anno_dir, "HMBA_BG_MWB_homology.xlsx"), sheet="MWB_Consensus_homology")
rownames(homology) = homology$Group
homology = homology[cluster_meta$Group, ]

## Fix cluster level homology to look good in sankey
homology$curated_ABC_WMB_cluster[is.na(homology$curated_ABC_WMB_cluster)] = homology$curated_ABC_WMB_supertype[is.na(homology$curated_ABC_WMB_cluster)]

homology_to_viz = homology %>% filter(!(Group %in% c("OPC", "COP", "lmOligo", "Oligo PLEKHG1", "Oligo OPALIN", "lmAstro",
                                                    "Astrocyte", "Ependymal", "Microglia", "BAM", "Monocyte", "B cells", "T cells", "VLMC",
                                                    "SMC", "Pericyte", "Endo", "SN SEMA5A GABA", "SN EBF2 GABA", "ZI-HTH GABA")))

##
homology_plot = homology_to_viz[,c(1,6,5,4)]

##
homology_plot$color = homology_plot$Group
homology_plot <- to_lodes_form(homology_plot, 
                                    value="stratum", 
                                    id="alluvium",
                                    axes = c(1:4))
str(homology_plot)
homology_plot$stratum = factor(homology_plot$stratum, levels=c(homology %>% arrange(factor(homology$Group, levels=cluster_meta$Group)) %>% distinct(Group) %>% pull(Group), 
                                                                Reduce(union, list(#homology %>% arrange(Group, curated_ABC_WMB_cluster) %>% distinct(curated_ABC_WMB_cluster) %>% pull(curated_ABC_WMB_cluster), 
                                                                                    homology %>% arrange(factor(Group, levels=cluster_meta$Group), factor(curated_ABC_WMB_supertype, levels=unique(homology$curated_ABC_WMB_supertype))) %>% distinct(curated_ABC_WMB_supertype) %>% pull(curated_ABC_WMB_supertype),
                                                                                    homology %>% arrange(factor(Group, levels=cluster_meta$Group), factor(curated_ABC_WMB_supertype, levels=unique(homology$curated_ABC_WMB_supertype))) %>% distinct(curated_ABC_WMB_subclass) %>% pull(curated_ABC_WMB_subclass),
                                                                                    homology %>% arrange(factor(Group, levels=cluster_meta$Group), factor(curated_ABC_WMB_supertype, levels=unique(homology$curated_ABC_WMB_supertype)), curated_ABC_WMB_subclass) %>% distinct(curated_ABC_WMB_class) %>% pull(curated_ABC_WMB_class)))
                                                                ))

##
ggplot(homology_plot,
       aes(x = x, stratum = stratum, alluvium=alluvium)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow(aes(color=color, fill=color), lwd=1.0) +
    geom_stratum(aes(color=color), alpha = 1.0) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
    scale_fill_manual(
        values = setNames(cluster_meta$color_hex_group, cluster_meta$Group),
        breaks = cluster_meta$Group
    ) +
    scale_color_manual(
        values = setNames(cluster_meta$color_hex_group, cluster_meta$Group),
        breaks = cluster_meta$Group
    ) +
    theme_void() +
    theme(legend.position = "none")

##
ggsave(file.path(analysis_dir, "primate_mouse_bg_homology.pdf"), width=12, height=10, units="in", dpi=600)
