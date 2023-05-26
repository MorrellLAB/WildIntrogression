#!/usr/bin/env Rscript

library(ggplot2)

wbdc_obs_het_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_vcf/compare_heterozygosity/wbdc_bopa_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt"
wbdc_obs_het_imputed_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_vcf/compare_heterozygosity/wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.observed_heterozygosity.txt"

dom_obs_het_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_vcf/compare_heterozygosity/domesticated_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt"
dom_obs_het_imputed_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_vcf/compare_heterozygosity/domesticated_snps.polymorphic.filt_miss_het.phased.imputed.observed_heterozygosity.txt"

#------------
wbdc_obs_het_df <- read.delim(wbdc_obs_het_fp, header=TRUE)
wbdc_obs_het_df["dataset"] <- "filtered"
wbdc_obs_het_imputed_df <- read.delim(wbdc_obs_het_imputed_fp, header=TRUE)
wbdc_obs_het_imputed_df["dataset"] <- "imputed"

cat_wbdc_obs_het <- rbind(wbdc_obs_het_df, wbdc_obs_het_imputed_df)

ggplot(wbdc_obs_het_df, aes(x=X..3.sample, y=observed_heterozygosity)) +
  geom_point(colour="dodgerblue", alpha=0.6) +
  theme_classic()

ggplot(wbdc_obs_het_imputed_df, aes(x=X..3.sample, y=observed_heterozygosity)) +
  geom_point(colour="dodgerblue", alpha=0.6) +
  theme_classic()

ggplot(cat_wbdc_obs_het, aes(x=dataset, y=observed_heterozygosity)) +
  geom_boxplot() +
  theme_bw()

summary(wbdc_obs_het_df$observed_heterozygosity)
summary(wbdc_obs_het_imputed_df$observed_heterozygosity)

#----------------
dom_obs_het_df <- read.delim(dom_obs_het_fp, header=TRUE)
dom_obs_het_df["dataset"] <- "filtered"
dom_obs_het_imputed_df <- read.delim(dom_obs_het_imputed_fp, header=TRUE)
dom_obs_het_imputed_df["dataset"] <- "imputed"

cat_dom_obs_het <- rbind(dom_obs_het_df, dom_obs_het_imputed_df)

ggplot(dom_obs_het_df, aes(x=X..3.sample, y=observed_heterozygosity)) +
  geom_point(colour="dodgerblue", alpha=0.6) +
  theme_classic()

ggplot(dom_obs_het_imputed_df, aes(x=X..3.sample, y=observed_heterozygosity)) +
  geom_point(colour="dodgerblue", alpha=0.6) +
  theme_classic()

ggplot(cat_dom_obs_het, aes(x=dataset, y=observed_heterozygosity)) +
  geom_boxplot() +
  theme_bw()

summary(dom_obs_het_df$observed_heterozygosity)
summary(dom_obs_het_imputed_df$observed_heterozygosity)
