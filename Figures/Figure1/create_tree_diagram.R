library(data.tree)
library(openxlsx)
library(ggraph)
library(igraph)
library(dplyr)
library(anndata)

## Helpful locations which are assumed to already exist
annodir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
figure_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/Figures/Figure1"

##
consensus_anno = read.xlsx(file.path(annodir, "HMBA_BG_consensus_annotation.xlsx"), sheet="consensus_anno_pre-print")

##
consensus_anno = consensus_anno %>%
                  distinct(Neighborhood, Class, Subclass, Group, .keep_all = TRUE) %>%
                  select(Neighborhood, Class, Subclass, Group, display_order_group)
consensus_anno = consensus_anno[consensus_anno$display_order_group,]
# consensus_anno$Class = meta.data[consensus_anno$Group,]$Class
# consensus_anno$Subclass = meta.data[consensus_anno$Group,]$Subclass

## --------------------------------------
consensus_anno$pathString <- paste("Brain", 
                                   consensus_anno$Neighborhood,
                                   consensus_anno$Class,
                                   consensus_anno$Subclass, 
                                   consensus_anno$Group,
                                   sep = "/")
anno = as.Node(consensus_anno)

##
anno_graph = ToDataFrameNetwork(anno)

## 
plot.me = ggraph(anno_graph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  geom_node_text(aes(label=name, filter=leaf) , angle=90 , hjust=1, nudge_y = -0.05) +
  geom_node_point(alpha=0.6) +
  ylim(-.5, NA) +
  theme_void()
ggsave(paste0(figure_dir, "/hierarchy_tree.pdf"), width=12, height=5)

# ##
# leaves_order = anno$Get(function(node) node$isLeaf)
# leaves_order = names(leaves_order)[leaves_order]

# ##
# leaves_order = gsub("SM", "Pericytes/SM", leaves_order)
# leaves_order = gsub("Cap Endothelial", "Ven/Cap Endothelial", leaves_order)

# ## 
# rownames(consensus_anno) = consensus_anno$Group
# consensus_anno = consensus_anno[leaves_order,]
# consensus_anno$hierarchy_diagram_order = 1:nrow(consensus_anno)
# consensus_anno = consensus_anno[,c("Class", "Subclass", "Group", "order", "hierarchy_diagram_order")]
# write.csv(consensus_anno, file=paste0(figure_dir, "consenus_spc_metadata.csv"))