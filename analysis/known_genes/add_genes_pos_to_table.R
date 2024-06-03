#!/usr/bin/env Rscript

library(dplyr)
library(readxl)
library(stringr)

setwd("~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions")

fp_table <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table_S2_jdf4_raw_merge_genbank_ids-clean.xlsx"
fp_dom_genes_bed <- "~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01/blast_out/filtered/dom_genes_combined.gene_seq.bed"

#-------------------
df_table <- read_excel(fp_table)
df_dom_genes <- read.delim(fp_dom_genes_bed, header=F, sep="\t")
colnames(df_dom_genes) <- c("Chrom", "Start", "End", "GeneID")
df_dom_genes$GenBankID <- str_split_fixed(df_dom_genes$GeneID, '_', 2)[,1]

df_combined <- full_join(df_table, df_dom_genes, by="GenBankID")

write_csv(df_combined, file="Table_S2_jdf4_raw_merge_genbank_ids.gene_regions.csv", na="")
