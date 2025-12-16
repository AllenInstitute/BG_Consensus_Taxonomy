## Load scrattch.mapping
library(anndata)
library(ggplot2)
library(ggrepel)
library(scrattch.hicat)

## Workind dir
setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia")

## Load in data
aibs.alignment = read.csv("Macaque_basalganglia_AIBS_mean_expression.csv")
nemo.alignment = read.csv("Macaque_basalganglia_NeMO_mean_expression.csv")

############################################################
## Data wrangling
############################################################

## Only compare genes we can match
common.genes = intersect(aibs.alignment[,1], nemo.alignment[,1])

############################################################
## Compare mean gene expression
############################################################

##
plot.me = data.frame(gene = common.genes, 
                        aibs=aibs.alignment[aibs.alignment[,1] %in% common.genes,2], 
                        optimus=nemo.alignment[nemo.alignment[,1] %in% common.genes,2],
                        diff = abs(aibs.alignment[aibs.alignment[,1] %in% common.genes,2]-nemo.alignment[nemo.alignment[,1] %in% common.genes,2]))
rownames(plot.me) = common.genes

##
plot.gg = ggplot(plot.me, aes(x=aibs, y=optimus)) + 
            geom_point(aes(colour=diff, fill=diff)) + 
            geom_abline(slope=1, intercept=0) + 
            xlim(0, max(c(max(aibs.alignment[aibs.alignment[,1] %in% common.genes,2]), max(nemo.alignment[nemo.alignment[,1] %in% common.genes,2])))) +
            ylim(0, max(c(max(aibs.alignment[aibs.alignment[,1] %in% common.genes,2]), max(nemo.alignment[nemo.alignment[,1] %in% common.genes,2])))) +
            geom_label_repel(data = plot.me[plot.me$diff > 1,], aes(label=gene),
                box.padding   = 0.35, 
                point.padding = 0.5,
                segment.color = 'grey50')
ggsave(plot.gg, file="~/aibs_opt_gene_means.png", width=7, height=6, units="in", dpi=600)

##
write.table(plot.me[plot.me$diff > 1,], sep="\t", file="~/aibs_opt_gene_means.tsv", row.names=F)

############################################################
## Compare number of features detected per cell
############################################################

## Number of genes expressed per cell
aibs.count.gene = colSums(aibs.counts > 1)
optimus.count.gene = colSums(optimus.counts > 1)

##
plot.me = data.frame(aibs = aibs.count.gene, 
                        optimus = optimus.count.gene, 
                        diff = abs(aibs.count.gene - optimus.count.gene))

plot.gg = ggplot(plot.me, aes(x=aibs, y=optimus)) + 
            geom_point(aes(colour=diff, fill=diff)) + 
            geom_abline(slope=1, intercept=0) + 
            xlim(0, max(c(max(aibs.count.gene), max(optimus.count.gene)))) +
            ylim(0, max(c(max(aibs.count.gene), max(optimus.count.gene))))
ggsave(plot.gg, file="aibs_opt_gene_counts.png", width=7, height=6, units="in", dpi=600)