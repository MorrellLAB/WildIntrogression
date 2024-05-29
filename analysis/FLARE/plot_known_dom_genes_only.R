library(data.table)
library(readxl)
library(chromPlot)
library(dplyr)
library(tidyr)

all_genes_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes/domestication-related_genes.morex_v3.clean.csv"
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
chr_lengths_fp <- "~/GitHub/WildIntrogression/analysis/FLARE/morex_v3_chr_lengths.txt"

#snps_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.bed"
snps_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.vcf.gz"

#-----------------
# Centromere
centromere <- fread(centromere_fp)
colnames(centromere) <- c("Chrom", "Start")
centromere$End <- centromere$Start
centromere$Name <- "centromere"

# Chromosome lengths
tmp_chr_lengths <- fread(chr_lengths_fp)
chr_lengths <- data.frame(Chrom=tmp_chr_lengths$V1, Start=0, End=tmp_chr_lengths$V2, Name="contig")
chr_centromere <- rbind(centromere, chr_lengths) %>% arrange(Chrom, Start)

# Load SNPs
snps <- fread(snps_fp, skip="#CHROM")
snps_df <- data.frame(Chrom=snps$`#CHROM`, Start=snps$POS-1, End=snps$POS, ID=snps$ID)

# Gene CSV
all_genes_df <- read.csv(all_genes_fp, header=TRUE)

# chromPlot of known domestication related genes only
# color_snps <- "#33415c"
# color_genes <- "#343a40"

#snps$Colors <- color_snps
#snps <- snps %>% as.data.frame()

# Prepare df for plotting
# Noticed a few duplicates where abbreviation was flipped
# Remove these rows
p.all_genes_df <- all_genes_df %>%
  dplyr::select("Chrom", "Start", "End", "Gene_Abbr") %>%
  drop_na() %>%
  as.data.frame()
# Rename columns to work with chromPlot
colnames(p.all_genes_df) <- c("Chrom", "Start", "End", "Name")
# Add band color
#p.all_genes_df$Colors <- color_genes
pb.all_genes_df <- p.all_genes_df
colnames(pb.all_genes_df) <- c("Chrom", "Start", "End", "ID")

chromPlot(gaps=chr_centromere, bands=snps_df)

# Plot
chromPlot(gaps=chr_centromere, bands=pb.all_genes_df, figCols=7, title="All known domestication-related genes",
          stat=p.all_genes_df, statCol="Value", statName="Value", noHist=TRUE)

chromPlot(gaps=chr_centromere, bands=p.all_genes_df, figCols=7, title="All known domestication-related genes",
          stat=p.all_genes_df, statCol="Value", statName="Value", noHist=TRUE)

chromPlot(gaps=chr_centromere, bands=p.all_genes_df, figCols=7, title="All known domestication-related genes",
          stat=p.all_genes_df, statCol="Value", statName="Value", noHist=TRUE)
chromPlot(gaps=chr_centromere, bands=snps, figCols=7, title=samp)
