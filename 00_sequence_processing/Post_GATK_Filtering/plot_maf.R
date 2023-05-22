#!/usr/bin/env Rscript

library(ggplot2)

fp <- "~/Downloads/temp_msi/introgression_project/vcf_summary/dom_and_wild_snps_biallelic.callable.cap50x.maf"

df <- read.delim(fp, header=TRUE)

ggplot(df, aes(MAF)) +
  geom_density(fill="dodgerblue", colour="black", alpha=0.3) +
  theme_light()

ggplot(df, aes(MAF)) +
  geom_histogram(fill="dodgerblue", colour="black", alpha=0.3) +
  theme_light()
