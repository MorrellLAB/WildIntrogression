#!/usr/bin/env Rscript

library(dplyr)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/introgressed_overlap_dom_genes/")

# Table S2 containing positions of 50k SNPs and known domestication gene positions when available
fp_table_s2 <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table_S2_jdf4_raw_merge_genbank_ids.gene_regions.new_snp_pos.csv"

#------------------
# Remove disease resistance genes
table_s2 <- read.csv(fp_table_s2) %>% filter(is.na(Resistance_Gene))

table_s2_wGenes <- table_s2 %>% filter(Chrom != "NA")

table_s2_without_genes <- table_s2 %>% filter(is.na(Chrom))

# Create BED file for gene regions
genes_bed <- data.frame(chrom=table_s2_wGenes$Chrom, start=table_s2_wGenes$Start, end=table_s2_wGenes$End, id=table_s2_wGenes$GenBankID)
# Create BED file for 50k SNPs where no gene sequences are available
snps_bed <- data.frame(chrom=table_s2_without_genes$chr_new, start=table_s2_without_genes$pos_bp_new-1, end=table_s2_without_genes$pos_bp_new, id=table_s2_without_genes$SNP.molecular.marker)

df_genes_snps <- rbind(genes_bed, snps_bed) %>% arrange(chrom, start)

# Save to file
write.table(df_genes_snps, file="dom_genes_50k_snps_noDiseaseRes.bed", row.names=F, col.names=F, sep="\t", quote=F)
