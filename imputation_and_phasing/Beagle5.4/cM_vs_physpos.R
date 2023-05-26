#!/usr/bin/env Rscript

library(ggplot2)

plink_map_wbdc_geno_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_vcf/wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map"

df_pm_wbdc_geno <- read.delim(plink_map_wbdc_geno_fp, header=FALSE)

ggplot(df_pm_wbdc_geno, aes(x=V4, y=V3)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~V1)
