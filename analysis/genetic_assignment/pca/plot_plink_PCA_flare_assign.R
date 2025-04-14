#!/usr/bin/env Rscript
library(stringr)
library(ggplot2)
library(dplyr)
library(plotly)
library(ggfortify)
library(ggrepel)
library(readxl)
library(vegan)

# Genotype data NSGC and wild
fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca/dom_and_wbdc_geno_samples.pca.eigenvec"

# Genotype data wild only
fp_wild <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca/wbdc_geno_samples.pca.eigenvec"
fp_wild_percent <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca/wbdc_geno_samples.pca.eigenval"

# GBS data wild only
fp_gbs_wild <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca/wbdc_gbs.pca.eigenvec"
fp_gbs_wild_percent <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca/wbdc_gbs.pca.eigenval"

# Final list of wild introgressed samples in Table S1 after FLARE, ABBABABA, etc. analyses
introgressed_fp <- "~/GitHub/WildIntrogression/analysis/FLARE/sample_list_wild_introgressed_tableS1.txt"
out_dir <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment"

# File containing longitude and latitude
geog_coord_fp <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v4/Tables_and_Files_v3/Table S1.xlsx"

#--------------------
# Read in files
df <- read.delim(file = fp, header = F, sep = " ")
df_wild <- read.delim(file = fp_wild,  header = F, sep = " ")
wild_eigenval <- read.delim(file = fp_wild_percent, header = F)
df_introgressed <- read.delim(file = introgressed_fp, header = F)
df_geog_coord <- read_excel(geog_coord_fp, sheet = "Sheet1", skip = 1)

# GBS wild
df_gbs <- read.delim(file=fp_gbs_wild, header=F, sep=" ")
# Modify WBDC naming to match that in sample list
# GBS data samples are formatted as "WBDC-001"
# Sample list is formatted as "WBDC001"
df_gbs$V1 <- sub("-", "", df_gbs$V1)
df_gbs$V2 <- sub("-", "", df_gbs$V2)
gbs_wild_eigenval <- read.delim(file=fp_gbs_wild_percent, header=F)

setwd(out_dir)

#--------------------
# Genotype data NSGC and wild
# Domesticated vs wild vs wild_introgressed samples
# Add groups to use for coloring points
# First add placeholder column
df['group'] <- NA

for (i in 1:nrow(df)) {
  if (str_detect(string=df[i,1], pattern="CIho")) {
    # Group is domesticated
    df[i, "group"] <- "Domesticated"
  } else if (str_detect(string=df[i,1], pattern="PI")) {
    # Group is domesticated
    df[i, "group"] <- "Domesticated"
  } else if (str_detect(string=df[i,1], pattern="OWB")) {
    # Group is domesticated
    df[i, "group"] <- "Domesticated"
  } else if (str_detect(string=df[i,1], pattern="WBDC")) {
    # Group is wild
    if (df[i,1] %in% df_introgressed$V1) {
      # Wild sample is introgressed
      df[i, "group"] <- "Wild Introgressed"
    } else {
      # Group is wild
      df[i, "group"] <- "Wild"
    }
  } else {
    print("Prefix doesn't match for sample:")
    df[i, 1]
  }
}

# Plot PC1 vs PC2
ggplot(data=df %>% arrange(group), aes(V3, V4, colour=as.factor(group))) +
  geom_point(alpha=0.6) +
  geom_text_repel(aes(label=V1), data=subset(df, group == "Wild Introgressed"), max.overlaps= 20, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#dce0d9", "#4ea5ff", "#ff0000")) +
  xlab("PC1") + ylab("PC2") +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc2-dom_and_wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC3 vs PC2
ggplot(data=df %>% arrange(group), aes(V5, V4, colour=as.factor(group))) +
  geom_point(alpha=0.6) +
  geom_text_repel(aes(label=V1), data=subset(df, group == "Wild Introgressed"), max.overlaps= 20, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#dce0d9", "#4ea5ff", "#ff0000")) +
  xlab("PC3") + ylab("PC2") +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc3_vs_pc2-dom_and_wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC1 vs PC3
ggplot(data=df %>% arrange(group), aes(V3, V5, colour=as.factor(group))) +
  geom_point(alpha=0.6) +
  geom_text_repel(aes(label=V1), data=subset(df, group == "Wild Introgressed"), max.overlaps= 20, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#dce0d9", "#4ea5ff", "#ff0000")) +
  xlab("PC1") + ylab("PC3") +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc3-dom_and_wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

#--------------------
# Genotype data
# Wild samples only
# Add groups to use for coloring points
# First add placeholder column
df_wild['group'] <- NA

for (i in 1:nrow(df_wild)) {
  if (df_wild[i,1] %in% df_introgressed$V1) {
    # Group is introgressed
    df_wild[i, "group"] <- "Wild Introgressed"
  } else if (str_detect(string=df_wild[i,1], pattern="WBDC")) {
    # Group is wild
    df_wild[i, "group"] <- "Wild"
  } else {
    print("Prefix doesn't match for sample:")
    df_wild[i, 1]
  }
}

# axis labels
wpc1lab <- paste0("PC1 (", round(wild_eigenval$V1[1], digits=2), "%)")
wpc2lab <- paste0("PC2 (", round(wild_eigenval$V1[2], digits=2), "%)")
wpc3lab <- paste0("PC3 (", round(wild_eigenval$V1[3], digits=2), "%)")

# Plot PC1 vs PC2
ggplot(data=df_wild %>% arrange(group), aes(V3, V4, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_wild, group == "Wild Introgressed"), max.overlaps= 10, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(wpc1lab) + ylab(wpc2lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc2-wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC3 vs PC2
ggplot(data=df_wild %>% arrange(group), aes(V5, V4, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_wild, group == "Wild Introgressed"), max.overlaps= 10, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(wpc3lab) + ylab(wpc2lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc3_vs_pc2-wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC1 vs PC3
ggplot(data=df_wild %>% arrange(group), aes(V3, V5, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_wild, group == "Wild Introgressed"), max.overlaps= 10, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(wpc1lab) + ylab(wpc3lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc3-wbdc_geno_tableS1.jpg", width=10, height=6, units="in", dpi=300)

#--------------------
# GBS data
# Wild samples only
# Add groups to use for coloring points
# First add placeholder column
df_gbs['group'] <- NA

for (i in 1:nrow(df_gbs)) {
  if (df_gbs[i,1] %in% df_introgressed$V1) {
    # Group is introgressed
    df_gbs[i, "group"] <- "Wild Introgressed"
  } else if (str_detect(string=df_gbs[i,1], pattern="WBDC")) {
    # Group is wild
    df_gbs[i, "group"] <- "Wild"
  } else {
    print("Prefix doesn't match for sample:")
    df_gbs[i, 1]
  }
}

# axis labels
gbswpc1lab <- paste0("PC1 (", round(gbs_wild_eigenval$V1[1], digits=2), "%)")
gbswpc2lab <- paste0("PC2 (", round(gbs_wild_eigenval$V1[2], digits=2), "%)")
gbswpc3lab <- paste0("PC3 (", round(gbs_wild_eigenval$V1[3], digits=2), "%)")

# Plot PC1 vs PC2
ggplot(data=df_gbs %>% arrange(group), aes(V3, V4, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_gbs, group == "Wild Introgressed"), max.overlaps= 11, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(gbswpc1lab) + ylab(gbswpc2lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc2-wbdc_GBS_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC3 vs PC2
ggplot(data=df_gbs %>% arrange(group), aes(V5, V4, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_gbs, group == "Wild Introgressed"), max.overlaps= 11, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(gbswpc3lab) + ylab(gbswpc2lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc3_vs_pc2-wbdc_GBS_tableS1.jpg", width=10, height=6, units="in", dpi=300)

# Plot PC1 vs PC3
ggplot(data=df_gbs %>% arrange(group), aes(V3, V5, colour=as.factor(group))) +
  geom_point(alpha=0.7) +
  geom_text_repel(aes(label=V1), data=subset(df_gbs, group == "Wild Introgressed"), max.overlaps= 11, show.legend=F) +
  theme_bw() +
  scale_color_manual(values=c("#4ea5ff", "#ff0000")) +
  xlab(gbswpc1lab) + ylab(gbswpc3lab) +
  theme(legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
ggsave("plink_pca_pc1_vs_pc3-wbdc_GBS_tableS1.jpg", width=10, height=6, units="in", dpi=300)
