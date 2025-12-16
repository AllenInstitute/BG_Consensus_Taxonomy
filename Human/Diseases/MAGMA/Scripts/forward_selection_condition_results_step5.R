library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
filename_cond <- args[1]
filename_marg <- args[2]
pheno <- args[3]
ancestry <- args[4]
folder <- args[5]

## READ MARGINAL AND CONDITIONAL RESULTS:
df_cond <- read_table(filename_cond, skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_cond = P)
df_marg <- read_table(filename_marg, skip=4) %>% 
  subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>% 
  dplyr::rename(P_marg = P) %>% 
  arrange(P_marg)
n_sig <- ceiling(sqrt(nrow(df_cond)))
n_clusters <- nrow(df_marg)


## Combine the marginal and conditional results: 
df_comb <- left_join(df_cond, df_marg, by = "VARIABLE")
df_comb <- df_comb %>% 
  mutate(VarCode = rep(c("a", "b"),times=nrow(df_comb)/2)) %>%
  pivot_wider(names_from = VarCode, values_from=c(VARIABLE, P_cond, P_marg)) %>%
  mutate(PS_a=log10(P_cond_a)/log10(P_marg_a), PS_b=log10(P_cond_b)/log10(P_marg_b))
# assert order
stopifnot(df_comb$P_marg_a <= df_comb$P_marg_b)

select_list <- df_marg$VARIABLE[1:n_sig]
reverse_list <- df_comb$VARIABLE_a[df_comb$PS_a < 0.2 & df_comb$PS_b >= 0.2]
select_list <- select_list[select_list %in% reverse_list == F][-1]

## Conduct forward selection:
# The independent list starts with the most marginally significant cell type
indep_list <- c(toString(df_marg[1,1])) 
for (cur_var in select_list) {
  for (indep_var in indep_list) {
    # drop if the condition indicates so with any cell types in indep_list
    df_cur = df_comb[df_comb$VARIABLE_a==indep_var & df_comb$VARIABLE_b==cur_var,]
    if ((df_cur$PS_a >= 0.8 & df_cur$PS_b >= 0.8) | 
        (df_cur$PS_a >= 0.5 & df_cur$PS_b >= 0.5 & df_cur$P_cond_b < 0.05)) {
      if (tail(indep_list,1) == indep_var) indep_list <- c(indep_list, cur_var)
    } else break
  }
}
cat(indep_list)

## GET THE EQUIVALENCE WITH THE CLUSTER LABEL:
indep_list_clean <- gsub("Cluster_", "", indep_list)
indep_list_ids <- as.integer(indep_list_clean)


output_path <- paste0(
  ".../", pheno, "_automated_new/cluster_ids_",
  ancestry, "_", pheno, ".csv"
)
write.csv(data.frame(indep_list_ids), output_path, row.names = FALSE)
