#!/usr/bin/env Rscript

fam_fp <- "~/Downloads/temp_msi/temp_admixture/pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.fam"
k2q_fp <- "~/Downloads/temp_msi/temp_admixture/pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.2.Q"
k3q_fp <- "~/Downloads/temp_msi/temp_admixture/pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.3.Q"
k7q_fp <- "~/Downloads/temp_msi/temp_admixture/pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.7.Q"
k20q_fp <- "~/Downloads/temp_msi/temp_admixture/pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.20.Q"

fam_df <- read.table(fam_fp)
df_k2q <- read.table(k2q_fp)
df_k3q <- read.table(k3q_fp)
df_k7q <- read.table(k7q_fp)
df_20q <- read.table(k20q_fp)

add_pop_group <- function(fam_df, df_k) {
  # Combine individual IDs with admixture results for sorting
  df_kq_combined <- data.frame(pop_group=fam_df$V1, df_k)s
  # Sort by IDs
  df_kq_csorted <- df_kq_combined[order(df_kq_combined$pop_group), ]
  # Remove IDs to make plotting easier
  df_kq_csp <- df_kq_csorted[, 2:ncol(df_kq_csorted)]
  return(df_kq_csp)
}

# df_k7q_combined <- data.frame(pop_group=fam_df$V1, df_k7q)
# df_k7q_csorted <- df_k7q_combined[order(df_k7q_combined$pop_group), ]
# df_k7q_csp <- df_k7q_csorted[, 2:ncol(df_k7q_csorted)]

df_k20_csp <- add_pop_group(fam_df, df_20q)

barplot(t(as.matrix(df_k2q)), col=rainbow(2), xlab="Individual #", ylab="Ancestry", border=NA)
barplot(t(as.matrix(df_k3q)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)
barplot(t(as.matrix(df_k7q_csp)), col=rainbow(7), xlab="Individual #", ylab="Ancestry", border=NA)

barplot(t(as.matrix(df_k20_csp)), col=rainbow(20), xlab="Individual #", ylab="Ancestry", border=NA)
