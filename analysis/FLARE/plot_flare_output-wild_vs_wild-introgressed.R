#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(dplyr)
library(readxl)
library(gdata)
library(R.utils)
#library(ggbeeswarm)
#library(ggridges)
#library(hrbrthemes)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots_wild_vs_wild-introgressed")

# Filepaths to FLARE output VCFs
# wild-introgressed
fp_wild_introgressed <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.vcf.gz"
fp_wild_rand1 <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_not_introgressed-rand1_48_samp_out/wild_rand1_48samp_geno.flare.out.anc.vcf.gz"
fp_wild_rand2 <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_not_introgressed-rand2_48_samp_out/wild_rand2_48samp_geno.flare.out.anc.vcf.gz"
fp_wild_rand3 <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_not_introgressed-rand3_48_samp_out/wild_rand3_48samp_geno.flare.out.anc.vcf.gz"

# Chromosome lengths and centromere
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
chr_lengths_fp <- "~/GitHub/WildIntrogression/analysis/FLARE/morex_v3_chr_lengths.txt"

#--------------------
# Centromere
centromere <- fread(centromere_fp)
colnames(centromere) <- c("Chrom", "Start")
centromere$End <- centromere$Start
centromere$Name <- "centromere"

# Chromosome lengths
tmp_chr_lengths <- fread(chr_lengths_fp)
chr_lengths <- data.frame(Chrom=tmp_chr_lengths$V1, Start=0, End=tmp_chr_lengths$V2, Name="contig")
chr_centromere <- rbind(centromere, chr_lengths) %>% arrange(Chrom, Start)

#---------------------
prep_flare_output <- function(flare_out_vcf_fp) {
  flare_df <- fread(flare_out_vcf_fp, skip="#CHROM")
  
  # Get a list of sample names only
  sample_names <- colnames(flare_df)[10:ncol(flare_df)]
  
  df <- data.frame(matrix(ncol=8, nrow=0))
  colnames(df) <- c("chr", "start", "end", "id", "sample", "gt", "an1", "an2")
  for (samp in sample_names) {
    print(samp)
    samp_col_idx <- which(colnames(flare_df) == samp)
    # Select relevant columns from df
    tmp_df <- flare_df[, c(1:3)] %>%
      mutate(start=POS-1, sample_name=samp)
    tmp_samp <- flare_df[[samp_col_idx]]
    # Reformat sample genotypes/annotations
    tmp_split <- data.frame(do.call('rbind', strsplit(as.character(tmp_samp), ':', fixed=TRUE)))
    colnames(tmp_split) <- c("GT", "AN1", "AN2")
    reformatted_df <- cbind(tmp_df, tmp_split) %>%
      relocate(start, .before=POS)
    colnames(reformatted_df) <- c("chr", "start", "end", "id", "sample", "gt", "an1", "an2")
    # Append to df
    df <- rbind(df, reformatted_df)
  }
  
  # Preset codes for populations, from header line in VCF output from FLARE
  # ##ANCESTRY=<cultivar=0,breeding=1,landrace=2,uncertain=3,genetic=4,wild=5>
  df_anc <- df %>%
    mutate(
      ancestry = case_when(
        gt == "0|0" ~ "hom_reference",
        gt != "0|0" & an1 == an2 & an1 == 5 ~ "wild",
        gt == "1|1" & an1 == an2 & an1 == 4 ~ "genetic",
        gt == "1|1" & an1 == an2 & an1 == 3 ~ "uncertain",
        gt == "1|1" & an1 == an2 & an1 == 2 ~ "landrace",
        gt == "1|1" & an1 == an2 & an1 == 1 ~ "breeding",
        gt == "1|1" & an1 == an2 & an1 == 0 ~ "cultivar",
        gt == "0|1" | gt == "1|0" ~ "het",
        .default = "other"
      )
    )
  return(df_anc)
}

# Read in FLARE output files and reformat for plotting
df_wild_introgressed <- prep_flare_output(fp_wild_introgressed)
df_wild_rand1 <- prep_flare_output(fp_wild_rand1)
df_wild_rand2 <- prep_flare_output(fp_wild_rand2)
df_wild_rand3 <- prep_flare_output(fp_wild_rand3)
unique(df_wild_introgressed$ancestry)
unique(df_wild_rand1$ancestry)
unique(df_wild_rand2$ancestry)
# wild_rand3 did not have any evidence of wild introgression
# Won't include in plots below
unique(df_wild_rand3$ancestry)

# Calculate average number of IBS tracts
# These are based on unique classes in the "ancestry" column
# Calculate proportion of SNPs assigned to each ancestry category for each sample
prop_wild_introgressed <- df_wild_introgressed %>%
  select(sample, ancestry) %>%
  group_by(sample, ancestry) %>%
  dplyr::summarise(count_anc=n()) %>%
  group_by(sample) %>%
  dplyr::mutate(prop = count_anc / sum(count_anc)) %>%
  as.data.frame()

prop_wild_rand1 <- df_wild_rand1 %>%
  select(sample, ancestry) %>%
  group_by(sample, ancestry) %>%
  dplyr::summarise(count_anc=n()) %>%
  group_by(sample) %>%
  dplyr::mutate(prop = count_anc / sum(count_anc))

prop_wild_rand2 <- df_wild_rand2 %>%
  select(sample, ancestry) %>%
  group_by(sample, ancestry) %>%
  dplyr::summarise(count_anc=n()) %>%
  group_by(sample) %>%
  dplyr::mutate(prop = count_anc / sum(count_anc))

prop_wild_rand3 <- df_wild_rand3 %>%
  select(sample, ancestry) %>%
  group_by(sample, ancestry) %>%
  dplyr::summarise(count_anc=n()) %>%
  group_by(sample) %>%
  dplyr::mutate(prop = count_anc / sum(count_anc))

# Prepare introgressed regions based on consecutive markers to save to file
prep_introg_regions <- function(df_anc) {
  df_introgressed_regions <- df_anc %>%
    group_by(rleid(ancestry), ancestry, chr, sample, gt, an1, an2) %>%
    summarise(start = gdata::first(start), end=gdata::last(end), .groups='drop') %>%
    dplyr::select(chr, start, end, sample, gt, an1, an2, ancestry) %>%
    as.data.frame()
  # Add size of introgressed intervals
  df_introgressed_regions$bp_length <- df_introgressed_regions$end - df_introgressed_regions$start
  df_introgressed_regions$Mbp_length <- df_introgressed_regions$bp_length / 1000000
  return(df_introgressed_regions)
}

regions_wild_introgressed <- prep_introg_regions(df_wild_introgressed)
regions_wild_rand1 <- prep_introg_regions(df_wild_rand1)
regions_wild_rand2 <- prep_introg_regions(df_wild_rand2)
regions_wild_rand3 <- prep_introg_regions(df_wild_rand3)
# Save to file
#write_delim(as.data.frame(df_introgressed_regions), file="wbdc_likely_introgressed_segments.txt", delim="\t", quote="none")

# Pull out categories for wild vs wild-introgressed IBS segment lengths comparison
wild_introgressed_breeding <- regions_wild_introgressed[regions_wild_introgressed$ancestry == "breeding", ]
wild_rand1_wild <- regions_wild_rand1[regions_wild_rand1$ancestry == "wild", ]
wild_rand2_wild <- regions_wild_rand2[regions_wild_rand2$ancestry == "wild", ]
wild_rand3_wild <- regions_wild_rand3[regions_wild_rand3$ancestry == "wild", ]

# Add column indicating grouping
wild_introgressed_breeding$Dataset <- "Wild-introgressed"
wild_rand1_wild$Dataset <- "Wild_Set1"
wild_rand2_wild$Dataset <- "Wild_Set2"
wild_rand3_wild$Dataset <- "Wild_Set3"

wild_introgressed_breeding$Type <- "Wild-introgressed"
wild_rand1_wild$Type <- "Wild"
wild_rand2_wild$Type <- "Wild"
wild_rand3_wild$Type <- "Wild"

# Combine into single df
df_comb <- rbind(wild_introgressed_breeding, wild_rand1_wild, wild_rand2_wild, wild_rand3_wild)

#-----------------
# Define colors for groups
color_breeding <- "#a70b0b"
color_wild <- "#bfdbf7"
#color_wild <- "#dee2e6"
color_hom_ref <- "#dee2e6"
color_genes <- "#353535"

#-----------------
# Plot distributions of IBS tract lengths
ggplot(df_comb, aes(x=Dataset, y=Mbp_length)) +
  geom_jitter(aes(colour=Type), alpha=0.6, width=0.35) +
  stat_summary(aes(group=Dataset), fun="mean", geom="crossbar", color="grey25", width=0.7, lwd=0.3) +
  stat_summary(aes(label=round(after_stat(y), 3)), fun="mean", colour="grey25", size=4, fontface="bold", geom="text", position=position_nudge(x=0.46)) +
  scale_color_manual(values=c(color_wild, color_breeding)) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  xlab("") +
  ylab("IBS Segment Length (Mbp)")

# Save plot
ggsave("ibs_lengths_wild_vs_wild-introgressed_jitterplot.jpg", width=12, height=10, units="in", dpi=300)

ggplot(df_comb, aes(x=chr, y=Mbp_length, fill=Dataset)) +
  geom_boxplot(alpha=0.8) +
  scale_fill_manual(values=c("#a2d6f9", "#1e96fc", "#072ac8", color_breeding)) +
  theme_bw() +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  xlab("") +
  ylab("IBS Segment Length (Mbp)")

# Save plot
ggsave("ibs_lengths_wild_vs_wild-introgressed_by_chr_boxplot.jpg", width=14, height=12, units="in", dpi=300)
 