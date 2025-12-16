library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
filename_marg_result_R <- args[1]
filename_output_R <- args[2]
genome_built <- args[3]
print("genome_built")
print(genome_built)
base_dir_BG_dataset <- args[4]

filename_spe_matrix = paste0(base_dir_BG_dataset, genome_built, "/conti_specificity_matrix.txt")

df <- read_tsv(filename_spe_matrix)
colnames(df)[1] <- "GENE"
print(head(df))
df_marg <- read_table(filename_marg_result_R, skip=4) %>% 
  arrange(P) %>%
  filter(P < 0.05/n())
  #filter(P < 0.0001)    # get all of the significant clusters after Bonferroni correction. 
print(head(df_marg))

sig_list <- df_marg$VARIABLE
print("sig_list:")
print(sig_list)
if (length(sig_list)==0){
  print("Not significant clusters")
}
df_sig <- df[c("GENE", sig_list)]
print(head(df_sig))
write_tsv(df_sig, filename_output_R)

