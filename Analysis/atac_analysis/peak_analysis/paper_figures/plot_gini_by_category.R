library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/ATAC-conservation/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/"

## spinalcord metadata
cluster_meta = read.xlsx("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables/consenus_spc_metadata.xlsx")
cluster_meta = cluster_meta[!duplicated(cluster_meta$Group), ]
cluster_meta$Group = gsub(" ", "_", cluster_meta$Group)

## Load combined counts
combined_counts = read.csv(file.path(data_dir, "peak_counts_by_celltype_and_species_biased_level3.csv"), row.names=1)

## Create a 3 column data with with cell type, species, and num peaks
biased_counts = reshape2::melt(combined_counts, id.vars=c("max_accessible_celltype"), variable.name="species", value.name="num_peaks")

biased_counts$max_accessible_celltype = factor(biased_counts$max_accessible_celltype, levels=gsub("/", "-", gsub(" ", "_", cluster_meta$Group)))
biased_counts$max_accessible_celltype

## Remove low n counts for species
biased_counts[which(biased_counts$species == "human" & biased_counts$max_accessible_celltype == "Visceral_Motor_Neurons"), "num_peaks"] = 0
biased_counts[which(biased_counts$species == "level3_num_peaks" & biased_counts$max_accessible_celltype == "Visceral_Motor_Neurons"), "num_peaks"] = 0

## Plot biased counts 
plot_data = biased_counts
plot_data$num_peaks = log(as.numeric(plot_data$num_peaks)+1, base=2)

plot_data[which(plot_data$species == "level3_num_peaks"), "num_peaks"] = -1 * plot_data[which(plot_data$species == "level3_num_peaks"), "num_peaks"]
ggplot(plot_data, aes(fill=species, y=num_peaks, x=max_accessible_celltype)) + 
            geom_bar(position="stack", stat="identity", color="black", linewidth=0.25) +
            scale_fill_manual(values=c("blue", "green", "orange", "grey"), breaks = c("human", "macaque", "mouse", "level3_num_peaks")) +
            theme(axis.text = element_text(color = "#000000"),
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  legend.position = "right",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.background = element_rect(fill = "#FFFFFF")) +
            ylim(-10,20)
ggsave(file.path(data_dir, "figures/peak_counts_by_celltype_and_species_biased_level3.pdf"), width=6, height=4, units="in")