#!/usr/bin/env Rscript
# Merge revisions of genes table

library(dplyr)
library(readxl)

setwd("~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions")

# Files to merge
fp_older_genbank_id <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table S2 jdf revised names 300524_ahs _jdf3_CL.xlsx"
# Revised jdf4 with typos fixed by CL
fp_revised_jdf4 <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table S2 jdf revised names 300524_ahs _jdf4_CL.xlsx"

#--------
df_genbank_id <- read_excel(fp_older_genbank_id, sheet="Table S2 jdf GenbankID merging")
df_genbank_id$`Position (bp)` <- as.character(df_genbank_id$`Position (bp)`)

df_revised <- read_excel(fp_revised_jdf4, sheet="TableS2_jdf_for_merge")

# Do a rough full merge, then cleanup in Excel after writing
df_merged <- full_join(df_genbank_id, df_revised, by=c("SNP molecular marker", "Chrom.", "Position (bp)", "Position (cM)", "Locus symbol")) %>% arrange(`No.`)

write_csv(df_merged, file="Table_S2_jdf4_raw_merge_genbank_ids.csv", na="")
