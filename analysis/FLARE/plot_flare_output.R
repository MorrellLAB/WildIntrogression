#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(dplyr)
library(readxl)
library(gdata)
library(R.utils)
library(chromPlot)
library(ggbeeswarm)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots")

flare_out_vcf_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.vcf.gz"
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
chr_lengths_fp <- "~/GitHub/WildIntrogression/analysis/FLARE/morex_v3_chr_lengths.txt"

# Spreadsheet containing known domestication related genes and bed positions
genes_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes/domestication-related_genes.morex_v3.csv"
# BED file of gene overlaps with introgressed regions with GenbankID as column 4
genes_intro_bed_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/known_dom_genes_overlap_introgressed_regions.bed"

#--------------------
flare_df <- fread(flare_out_vcf_fp, skip="#CHROM")

# Centromere
centromere <- fread(centromere_fp)
colnames(centromere) <- c("Chrom", "Start")
centromere$End <- centromere$Start
centromere$Name <- "centromere"

# Chromosome lengths
tmp_chr_lengths <- fread(chr_lengths_fp)
chr_lengths <- data.frame(Chrom=tmp_chr_lengths$V1, Start=0, End=tmp_chr_lengths$V2, Name="contig")
chr_centromere <- rbind(centromere, chr_lengths) %>% arrange(Chrom, Start)

# Load CSV with gene info
genes_df <- read.csv(genes_fp, header=TRUE)
# Load BED file with gene related info
bed_df <- read.delim(genes_intro_bed_fp, sep="\t", header=FALSE)
colnames(bed_df) <- c("Chrom", "Start", "End", "GenBank_ID")
# Split names, some are formatted as genbank_id-transcript_id
# Example: FJ974009.1-HORVU.MOREX.r3.7HG0637550.1
bed_df <- bed_df %>%
  separate(col=GenBank_ID, into=c("GenBank_ID", "Transcript_ID"), sep="-", fill="right") %>%
  dplyr::select(-c("Transcript_ID"))
# Add gene abbreviations to bed df
intro.genes <- left_join(bed_df, genes_df, by=join_by("GenBank_ID", "Chrom", "Start", "End"))
# Save to file
write_csv(intro.genes, file="known_dom_genes_introgressed_regions.csv", quote="none", col_names=TRUE)

#---------------------
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

#---------------------
# Plot proportion of SNPs assigned to each group
# Calculate proportion of SNPs assigned to each group
prop_anc_summary <- df_anc %>%
  group_by(sample) %>%
  count(ancestry) %>%
  mutate(anc_prop = n / sum(n))
# Re-order groups
unique(prop_anc_summary$ancestry)
prop_anc_summary$ancestry <- factor(prop_anc_summary$ancestry, levels=c("breeding", "wild", "hom_reference"))
# Save to file
write_delim(prop_anc_summary, file="wild_introgressed_prop.txt", delim="\t", quote="none", col_names=TRUE)

mean(as.data.frame(prop_anc_summary[prop_anc_summary$ancestry == "breeding", "anc_prop"])$anc_prop)
sd(as.data.frame(prop_anc_summary[prop_anc_summary$ancestry == "breeding", "anc_prop"])$anc_prop)

# Check no introgression individuals
prop_snps_breeding <- prop_anc_summary[prop_anc_summary$ancestry == "breeding", ] %>%
  arrange(anc_prop)
# Identify individuals with <1%
lt1percent <- prop_snps_breeding[prop_snps_breeding$anc_prop < 0.01, ] %>% arrange(sample)
# Save to file
write_delim(lt1percent, file="wild_introgressed_prop.lt1percent.txt", delim="\t", quote="none", col_names=TRUE)
# Identify individuals with <2%
lt2percent <- prop_snps_breeding[prop_snps_breeding$anc_prop < 0.02, ] %>% arrange(sample)
# Save to file
write_delim(lt1percent, file="wild_introgressed_prop.lt2percent.txt", delim="\t", quote="none", col_names=TRUE)

# Sorted by sample names
ggplot(prop_anc_summary, aes(x=sample, y=anc_prop, fill=ancestry)) +
  geom_col(position=position_dodge2(width = 0.7, preserve = "single")) +
  scale_fill_manual(values=c("#476C9B", "#ADD9F4", "#dee2e6")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Wild Samples") +
  ylab("Proportion of SNPs in Group")
# Save plot
ggsave("wild_introgressed_prop.jpg", width=18, height=6, units="in", dpi=300)

#-----------------
# Plot per sample introgressed segments on chromosomes
# Define colors for groups
color_breeding <- "#a70b0b"
color_wild <- "#bfdbf7"
#color_wild <- "#dee2e6"
color_hom_ref <- "#dee2e6"
color_genes <- "#353535"

# Prepare introgressed regions based on consecutive markers to save to file
df_introgressed_regions <- df_anc %>%
  group_by(rleid(ancestry), ancestry, chr, sample, gt, an1, an2) %>%
  summarise(start = gdata::first(start), end=gdata::last(end), .groups='drop') %>%
  dplyr::select(chr, start, end, sample, gt, an1, an2, ancestry) %>%
  as.data.frame()
# Add size of introgressed intervals
df_introgressed_regions$bp_length <- df_introgressed_regions$end - df_introgressed_regions$start
df_introgressed_regions$Mbp_length <- df_introgressed_regions$bp_length / 1000000
# Save to file
write_delim(as.data.frame(df_introgressed_regions), file="wbdc_likely_introgressed_segments.txt", delim="\t", quote="none")

#-----------------
# Calculate average number and size of introgressed regions
df_breeding_introgressed <- df_introgressed_regions[df_introgressed_regions$ancestry == "breeding", ]
# Save to file
write_csv(df_breeding_introgressed, file="wbdc_likely_introgressed_segments-breeding.csv", quote="none", col_names=TRUE)

# Proportion of genome by combined total amount introgressed
prop_total_introg <- df_breeding_introgressed %>%
  group_by(sample) %>%
  dplyr::summarise(count_introgressed=n(),
                   total_Mbp_length_sample=sum(Mbp_length),
                   prop_of_genome_per_sample=total_Mbp_length_sample/4200) %>%
  as.data.frame()

mean(prop_total_introg$prop_of_genome_per_sample)
sd(prop_total_introg$prop_of_genome_per_sample)
max(prop_total_introg$prop_of_genome_per_sample)

prop_total_introg %>%
  arrange(prop_of_genome_per_sample)

ggplot(prop_total_introg, aes(x=prop_of_genome_per_sample, y=count_introgressed)) +
  geom_point() +
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  xlab("Proportion of the genome introgressed\n(per sample)") +
  ylab("Number of introgressed segments")
# Save plot
ggsave("prop_introgresed_vs_num_introgressed.jpg", width=10, height=8, units="in", dpi=300)

# Check count of SNPs lt1%
lt1percent %>%
  mutate(num_snps_lt1p=n*anc_prop)
# Check size of segments in these samples
df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC066", ] %>% dplyr::summarise(mean(Mbp_length))
df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC125", ] %>% dplyr::summarise(mean(Mbp_length))
df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC128", ] %>% dplyr::summarise(mean(Mbp_length))
df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC231", ] %>% dplyr::summarise(mean(Mbp_length))
# lt1% samples
lt1p_samp <- rbind(df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC066", ],
      df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC125", ],
      df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC128", ],
      df_breeding_introgressed[df_breeding_introgressed$sample == "WBDC231", ])
lt1p_samp %>% dplyr::summarise(mean(Mbp_length))
lt1p_samp %>% dplyr::summarise(mean(bp_length))
# Avg bp length lt1p
2490426 / 516505932
2490426 / sum(chr_lengths$End)

# Total number of tracts
nrow(df_breeding_introgressed)

# Format regions that are "breeding" to BED
bed_introgressed <- data.frame(chr=df_breeding_introgressed$chr, start=df_breeding_introgressed$start-1, end=df_breeding_introgressed$end, sample=df_breeding_introgressed$sample)
# Save to file
write_delim(bed_introgressed, file="wbdc_likely_introgressed_segments-breeding.bed", delim="\t", quote="none", col_names=FALSE)

# Per individual counts
df_breeding_introgressed %>%
  group_by(sample) %>%
  dplyr::summarize(num_introgressed=n()) %>%
  write_csv(file="per_sample_counts_of_introgressed_regions.csv", col_names=TRUE, quote="none")

df_breeding_introgressed %>%
  group_by(sample) %>%
  dplyr::summarize(num_introgressed=n()) %>%
  dplyr::summarize(mean(num_introgressed))

# Average length
mean(df_breeding_introgressed$Mbp_length)
sd(df_breeding_introgressed$Mbp_length)

# Per chromosome average length
df_breeding_introgressed %>%
  group_by(chr) %>%
  dplyr::summarize(Mean=mean(Mbp_length), Standard_deviation=sd(Mbp_length)) %>%
  write_csv(file="mean_per_chr_introgressed_length_breeding.csv", col_names=TRUE, quote="none")

# Check size distribution
summary(df_breeding_introgressed$Mbp_length)
mbp_percentile <- quantile(df_breeding_introgressed$Mbp_length, probs=c(0, 0.01, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.996, 0.997, 0.998, 0.999, 1)) %>%
  as.data.frame()
colnames(mbp_percentile) <- "Mbp_length"

sd(df_breeding_introgressed$Mbp_length)
mean(df_breeding_introgressed$Mbp_length)
sd(df_breeding_introgressed$Mbp_length)*2

# 95th percentile
nrow(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 7.5, ])
# 35 samples
length(unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 7.5, "sample"]))
# 18 samples (99th percentile)
length(unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 21, "sample"]))
unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 21 & df_breeding_introgressed$Mbp_length <= 49, "sample"])
# 9 samples (99.6th percentile)
length(unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length > 49, "sample"]))
unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length > 49, "sample"])
# 6 samples
length(unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 100, "sample"]))

setdiff(unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 21 & df_breeding_introgressed$Mbp_length <= 49, "sample"]), unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length > 49, "sample"]))

# Number of segments > 50 Mbp
nrow(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 50, ])
nrow(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 100, ])
nrow(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 200, ])

df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 100, ]
df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 200, ]
unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 100, c("chr", "start", "end", "sample")])
unique(df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 50, c("chr", "start", "end", "sample")])

df_breeding_introgressed[df_breeding_introgressed$Mbp_length >= 50, ]

#-----------------

# Plot length distribution of wild-domesticated introgression tracts
ggplot(df_breeding_introgressed, aes(x=chr, y=Mbp_length)) +
  geom_boxplot(outlier.shape = 1) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  xlab("Chromosome") +
  ylab("Introgressed Segment Size (Mbp)")
# Save plot
ggsave("wild_introgressed_size_distribution.jpg", width=10, height=6, units="in", dpi=300)

# Prepare version for chromPlot
df_regions_gt <- df_anc %>%
  group_by(rleid(ancestry), ancestry, chr, sample, gt) %>%
  summarise(start = gdata::first(start), end=gdata::last(end), .groups='drop') %>%
  dplyr::select(chr, start, end, sample, gt, ancestry) %>%
  mutate(
    Colors = case_when(
      ancestry == "hom_reference" ~ color_hom_ref,
      ancestry == "wild" ~ color_wild,
      ancestry == "breeding" ~ color_breeding
    )
  ) %>%
  as.data.frame()

# Rename columns
colnames(df_regions_gt) <- c("Chrom", "Start", "End", "sample", "gt", "Name", "Colors")

# Prepare gene info for plotting
p.intro.genes <- intro.genes %>%
  dplyr::select(Chrom, Start, End, Gene_Abbr)
colnames(p.intro.genes) <- c("Chrom", "Start", "End", "ID")
# Add column of colors
p.intro.genes$Colors <- color_genes

setDT(p.intro.genes)

for (samp in sample_names) {
  # Select sample
  curr_sample <- df_regions_gt[df_regions_gt$sample == samp, ] %>%
    dplyr::select(Chrom, Start, End, Name, Colors)
  # Pull genes that overlap introgressed segments only
  tmp_dom_curr_samp <- curr_sample[curr_sample$Name == "breeding", ]
  setDT(tmp_dom_curr_samp)
  setkey(tmp_dom_curr_samp, Chrom, Start, End)
  # Get region overlaps
  tmp_intro_genes <- foverlaps(p.intro.genes, tmp_dom_curr_samp, type="any") %>%
    filter(Start != "NA") %>%
    dplyr::select(Chrom, i.Start, i.End, ID, i.Colors) %>%
    as.data.frame()
  colnames(tmp_intro_genes) <- c("Chrom", "Start", "End", "ID", "Colors")
  # Generate and save plot
  out_file <- paste0(samp, "_chrom_plot.jpg")
  jpeg(out_file, units="in", width=7, height=3, res=300)
  print(samp)
  if (nrow(tmp_intro_genes) > 1) {
    # There are genes that overlap introgressed regions for current sample
    chromPlot(gaps=chr_centromere, bands=curr_sample, figCols=7, title=samp,
              stat=tmp_intro_genes, statCol="Value", statName="Value", noHist=TRUE)
  } else {
    chromPlot(gaps=chr_centromere, bands=curr_sample, figCols=7, title=samp)
  }
  #legend("bottomright", legend=c("breeding", "wild", "hom_reference"),
  #       col=c(color_breeding, color_wild, color_hom_ref), lty=1, lwd=5, cex=0.9)
  dev.off()
}
