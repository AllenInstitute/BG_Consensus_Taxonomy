## Code is adapted from Bryois (2020)and Duncan et al. The script produces continuous specificity MAGMA inputs. ##


library(tidyverse)
library(rhdf5)
library(AnnotationDbi)
library(dplyr)
library(tidyr)

genome_built = "NCBI37"

base_dir_BG_dataset = '.../data/newest_BG/'
file = paste0(base_dir_BG_dataset,'newest_BG_log2.h5')

h5f <- H5Fopen(file)
exp <- as.data.frame(t(h5f$matrix))

# Add cluster names 
n_clusters <- ncol(exp)  
level_name <- "Cluster_" 
clusters <- paste0(level_name, 0:(n_clusters - 1))
colnames(exp) <- clusters

# Read gene names
gene_names <- h5read(h5f, "gene_ids")
gene_names <- as.character(gene_names)

# Set row names of expression matrix
rownames(exp) <- gene_names
exp$Gene <- rownames(exp)

h5closeAll()


exp <- exp %>% add_count(Gene) %>% 
  filter(n==1) %>%
  dplyr::select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()

## LOAD GENE COORDINATES:
gene_coordinates <- 
  read_tsv(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Alma/data/MAGMA_external/", genome_built, ".3.gene.loc.extendedMHCexcluded"),
           col_names = c("ENTREZ", "chr", "start", "end"),col_types = 'cciicc') %>%
  dplyr::select(1:4) %>%
  mutate(chr=paste0("chr",chr))


## TABLE FOR ENTREZ AND ENSEMBL gene names: 
entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
    # Only keep genes with a unique ENtrez and ensembl id: 
entrez_ensembl_unique_genes_entrez <- entrez_ensembl %>% count(gene_id) %>% filter(n==1)
entrez_ensembl_unique_genes_ens <- entrez_ensembl %>% count(ensembl_id) %>% filter(n==1)
entrez_ensembl <- filter(entrez_ensembl,gene_id%in%entrez_ensembl_unique_genes_entrez$gene_id & ensembl_id %in% entrez_ensembl_unique_genes_ens$ensembl_id)
colnames(entrez_ensembl)[1] <- "ENTREZ"
colnames(entrez_ensembl)[2] <- "Gene"
gene_coordinates <- inner_join(entrez_ensembl,gene_coordinates) %>% as.tibble()


## WHEN ensembl is already included:
exp_annotated <- inner_join(exp, gene_coordinates, by = "Gene")

##Remove genes not expressed
exp_agg <- exp_annotated %>% dplyr::rename(ClusterID=column, Expr_sum_mean=Expr)

not_expressed <- exp_agg %>% 
  group_by(Gene) %>% 
  summarise(total_sum=sum(Expr_sum_mean)) %>% 
  filter(total_sum==0) %>% 
  dplyr::select(Gene) %>% unique() 

exp_agg <- filter(exp_agg,!Gene%in%not_expressed$Gene)
##Each cell type is scaled to the same total number of (transformed) molecules.
exp_agg <- exp_agg %>% 
  group_by(ClusterID) %>% 
  mutate(Expr_sum_mean=Expr_sum_mean*1000/sum(Expr_sum_mean, na.rm = TRUE)) %>% 
  ungroup()

exp_agg <- exp_agg %>%
  group_by(Gene) %>%
  mutate(sum_expr = sum(Expr_sum_mean, na.rm = TRUE)) %>%
  ungroup()

## Calculate Specificty: 
exp_agg <- exp_agg %>% 
  group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean, na.rm = TRUE)) %>%
  ungroup()


##Write MAGMA input file
exp_conti_spe <- exp_agg %>% dplyr::select(ENTREZ, ClusterID, specificity) %>% spread(ClusterID, specificity)
#colnames(exp_conti_spe)[1] <- "GENE"

exp_conti_spe %>% write_tsv(paste0(base_dir_BG_dataset, genome_built, "/conti_specificity_matrix.txt"))
