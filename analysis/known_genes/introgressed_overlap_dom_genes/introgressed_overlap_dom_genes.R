#!/usr/bin/env Rscript

library(dplyr)
library(readxl)
library(forcats)
library(ggplot2)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/introgressed_overlap_dom_genes")

#fp_bed <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/introgressed_segments_overlap_domestication_genes.bed"
fp_bed <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/introgressed_overlap_dom_genes/introgressed_segments_overlap_dom_genes_50k_snps_noDiseaseRes.bed"

fp2_bed <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/introgressed_overlap_dom_genes/introgressed_segments_overlap_dom_genes_noDiseaseRes.bed"

#------------------
df <- read.delim(fp_bed, header=F)
colnames(df) <- c("chr", "start", "end", "sample")

num_regions_overlap_genes_df <- df %>%
  group_by(sample) %>%
  dplyr::summarise(num_regions_overlap_genes=n()) %>%
  arrange(num_regions_overlap_genes) %>%
  mutate(sample=fct_reorder(sample, num_regions_overlap_genes))

ggplot(num_regions_overlap_genes_df, aes(x=sample, y=num_regions_overlap_genes)) +
  geom_col() +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) +
  ylim(0, 20) +
  xlab("Wild Samples") +
  ylab("Number of introgressed \nregions overlapping genes")
# Save plot
ggsave("wild_introgressed_overlap_dom_genes_50k_snps_noDiseaseRes.jpg", width=14, height=8, units="in", dpi=300)

#--------------------
df2 <- read.delim(fp2_bed, header=F)
colnames(df2) <- c("chr", "start", "end", "sample")

num_regions_overlap_genes_df2 <- df2 %>%
  group_by(sample) %>%
  dplyr::summarise(num_regions_overlap_genes=n()) %>%
  arrange(num_regions_overlap_genes) %>%
  mutate(sample=fct_reorder(sample, num_regions_overlap_genes))

ggplot(num_regions_overlap_genes_df2, aes(x=sample, y=num_regions_overlap_genes)) +
  geom_col() +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Wild Samples") +
  ylab("Number of introgressed \nregions overlapping genes")
# Save plot
ggsave("wild_introgressed_overlap_dom_genes_noDiseaseRes.jpg", width=14, height=8, units="in", dpi=300)
