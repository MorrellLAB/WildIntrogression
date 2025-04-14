#!/usr/bin/env Rscript

library(dplyr)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare")

# User provided input arguments
phys_pos_fp = "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.chr_pos_only.txt"
cM_fp = "~/Dropbox/Projects/Wild_Introgression/Analyses/recombination/bopa_and_9k_genetic_pos.txt"

#-------------
df_phys_pos = read.delim(phys_pos_fp, sep="\t", header=T)
colnames(df_phys_pos) <- c("chrom", "pos", "ID")

df_cM = read.delim(cM_fp, sep="\t", header=F) %>%
    mutate(V2=paste0("chr", V2))
colnames(df_cM) <- c("ID", "chrom", "cM")

df_combined = left_join(df_phys_pos, df_cM)

write_delim(x=df_combined, file="wild_likely_introgressed_geno.flare.out.anc.chr_pos_cM.txt", delim="\t", quote="none")
