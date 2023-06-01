#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(readr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# User provided input arguments
# hap-ibd.out.ibd.gz file
fp <- args[1]
out_dir <- args[2]
# Output files and plots will have this prefix
out_prefix <- args[3]

#--------------
df <- fread(fp, header=FALSE)
colnames(df) <- c("samp1id", "samp1hapidx", "samp2id", "samp2hapidx", "chr", "start_bp", "end_bp", "cM_length_ibd")
# Calculate bp size
df$bp_length_ibd <- df$end_bp - df$start_bp
df$Mbp_length_ibd <- df$bp_length_ibd / 1000000

# Identify rows where both samp1id and samp2id are WBDC
tmp.wbdc.wbdc <- df %>% filter(grepl("^WBDC", samp1id) & grepl("^WBDC", samp2id))
# Only include rows where WBDC is present in exactly either samp1id or samp2id but NOT both columns
df.wbdc <- df %>% filter(grepl("^WBDC", samp1id) | grepl("^WBDC", samp2id))
# Remove rows where both samp1id and samp2id are WBDC
# We are interested in WBDC and domesticated accession pairs only
df.filt <- setdiff(df.wbdc, tmp.wbdc.wbdc)

setwd(out_dir)
# Save to file
write_delim(df.filt, file=paste0(out_prefix, ".wbdc-dom_pairs_only.txt"), delim="\t", quote="none")

# Plot distribution of sizes
ggplot(df.filt, aes(x=chr, y=Mbp_length_ibd)) +
  geom_boxplot() +
  theme_classic() +
  xlab("Chromosome") +
  ylab("IBD segment length (Mbp)") +
  ylim(c(0, max(df.filt$Mbp_length_ibd))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
# Save plot to file
ggsave(paste0(out_prefix, "-Mbp_ibd_length.jpg"), width=10, height=6, units="in", dpi=300)

ggplot(df.filt, aes(x=chr, y=cM_length_ibd)) +
  geom_boxplot() +
  theme_classic() +
  xlab("Chromosome") +
  ylab("IBD segment length (cM)") +
  ylim(c(0, max(df.filt$cM_length_ibd))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key.size=unit(10, "point"))
# Save plot to file
ggsave(paste0(out_prefix, "-cM_ibd_length.jpg"), width=10, height=6, units="in", dpi=300)
