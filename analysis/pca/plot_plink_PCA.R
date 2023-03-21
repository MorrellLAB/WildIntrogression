#!/usr/bin/env Rscript
library(stringr)
library(ggplot2)

fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/dom_and_wbdc_samples.pca.eigenvec"
fp_wild <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/wbdc_samples.pca.eigenvec"
fp_wild_percent <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/wbdc_samples.pca.eigenval"
introgressed_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/likely_introgressed_sample_names.txt"
out_dir <- "~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment"

# Read in files
df <- read.delim(file = fp, header = F, sep = " ")
df_wild <- read.delim(file = fp_wild,  header = F, sep = " ")
wild_eigenval <- read.delim(file = fp_wild_percent, header = F)
df_introgressed <- read.delim(file = introgressed_fp, header = F)

#--------------------
# Domesticated and wild samples
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
    df[i, "group"] <- "Wild"
  } else {
    print("Prefix doesn't match for sample:")
    df[i, 1]
  }
}

# Plot PC1 vs PC2
pc1_pc2_out_filename <- paste0(out_dir, "/plink_pca_pc1_vs_pc2.pdf")
ggplot(data=df, aes(V3, V4, colour=as.factor(group))) +
  geom_point(alpha=0.3) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  theme(legend.title=element_blank())
ggsave(pc1_pc2_out_filename)

# Plot PC3 vs PC2
pc3_pc2_out_filename <- paste0(out_dir, "/plink_pca_pc3_vs_pc2.pdf")
ggplot(data=df, aes(V5, V4, colour=as.factor(group))) +
  geom_point(alpha=0.3) +
  theme_bw() +
  xlab("PC3") + ylab("PC2") +
  theme(legend.title=element_blank())
ggsave(pc3_pc2_out_filename)

# Plot PC1 vs PC3
pc1_pc3_out_filename <- paste0(out_dir, "/plink_pca_pc1_vs_pc3.pdf")
ggplot(data=df, aes(V3, V5, colour=as.factor(group))) +
  geom_point(alpha=0.3) +
  theme_bw() +
  xlab("PC1") + ylab("PC3") +
  theme(legend.title=element_blank())
ggsave(pc1_pc3_out_filename)

#--------------------
# Wild samples only
# Add groups to use for coloring points
# First add placeholder column
df_wild['group'] <- NA

for (i in 1:nrow(df_wild)) {
  if (df_wild[i,1] %in% df_introgressed$V1) {
    # Group is introgressed
    df_wild[i, "group"] <- "Introgressed"
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
w_pc1_pc2_out_filename <- paste0(out_dir, "/plink_pca_wild_only_pc1_vs_pc2.pdf")
ggplot(data=df_wild, aes(V3, V4, colour=as.factor(group))) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab(wpc1lab) + ylab(wpc2lab) +
  theme(legend.title=element_blank())
ggsave(w_pc1_pc2_out_filename)

# Plot PC3 vs PC2
w_pc3_pc2_out_filename <- paste0(out_dir, "/plink_pca_wild_only_pc3_vs_pc2.pdf")
ggplot(data=df_wild, aes(V5, V4, colour=as.factor(group))) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab(wpc3lab) + ylab(wpc2lab) +
  theme(legend.title=element_blank())
ggsave(w_pc3_pc2_out_filename)

# Plot PC1 vs PC3
w_pc1_pc3_out_filename <- paste0(out_dir, "/plink_pca_wild_only_pc1_vs_pc3.pdf")
ggplot(data=df_wild, aes(V3, V5, colour=as.factor(group))) +
  geom_point(alpha=0.5) +
  theme_bw() +
  xlab(wpc1lab) + ylab(wpc3lab) +
  theme(legend.title=element_blank())
ggsave(w_pc1_pc3_out_filename)
