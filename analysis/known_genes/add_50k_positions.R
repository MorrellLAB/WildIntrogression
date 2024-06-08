#!/usr/bin/env Rscript

# Add our 50K SNP positions to Table S2

library(readxl)
library(dplyr)

setwd("~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions")

fp_our_50k_pos <- "~/GitHub/WildIntrogression/analysis/known_genes/table_S2_snps_50k_idt90_noRescuedSNPs.bed"
fp_table_s2 <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table S2 clean jdf positions.xlsx"

fp_table_s2_plotting <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table_S2_jdf4_raw_merge_genbank_ids.gene_regions.xlsx"

#-----------------
snp_50k <- read.delim(fp_our_50k_pos, sep="\t", header=FALSE) %>%
    dplyr::select(!(V2))

df_table_s2 <- read_excel(fp_table_s2, skip=1) %>%
    dplyr::filter(row_number() <= n()-7)

df_new_pos <- full_join(snp_50k, df_table_s2, by=c("V4" = "SNP molecular marker"))
df_new_pos$Chrom. <- paste("chr", df_new_pos$Chrom., sep="")

df_new_pos$chr_concordant <- df_new_pos$V1 == df_new_pos$Chrom.
df_new_pos$pos_concordant <- df_new_pos$V3 == df_new_pos$`Position (bp)`
df_new_pos$orig_minus_our_pos_bp <- as.numeric(df_new_pos$`Position (bp)`) - as.numeric(df_new_pos$V3)

df_new_pos_sorted <- df_new_pos %>% arrange(No., V1, V3)

colnames(df_new_pos_sorted)[1:3] <- c("chr_new", "pos_bp_new", "SNP molecular marker")

write_csv(df_new_pos_sorted, file="Table S2 new pos draft.csv", na="", col_names=TRUE)

# Add our 50k snp positions to spreadsheet for plotting
df_table_s2_plotting <- read_excel(fp_table_s2_plotting, sheet="Table_S2_annotations")

dfp_new_pos <- full_join(snp_50k, df_table_s2_plotting, by=c("V4" = "SNP molecular marker"))
colnames(dfp_new_pos)[1:3] <- c("chr_new", "pos_bp_new", "SNP molecular marker")

write_csv(dfp_new_pos, file="Table_S2_jdf4_raw_merge_genbank_ids.gene_regions.new_snp_pos.csv")
