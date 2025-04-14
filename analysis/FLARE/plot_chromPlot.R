#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(readxl)
library(gdata)
library(R.utils)
library(chromPlot)
library(stringr)

# Updated chromPlots only with revised Table S2 gene list
# setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots_updated_chromPlots_our_pos/")
setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots_updated_chromPlots_our_pos_bw/")

flare_out_vcf_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild_likely_introgressed_geno.flare.out.anc.vcf.gz"
centromere_fp <- "~/GitHub/morex_reference/morex_v3/pericentromere/MorexV3_centromere_positions.tsv"
chr_lengths_fp <- "~/GitHub/WildIntrogression/analysis/FLARE/morex_v3_chr_lengths.txt"

# Spreadsheet containing known domestication related genes and bed positions
genes_fp <- "~/Dropbox/Projects/Wild_Introgression/Jerry_gene_table_revisions/Table_S2_jdf4_raw_merge_genbank_ids.gene_regions.new_snp_pos.csv"
# BED file of gene overlaps with introgressed regions with GenbankID as column 4
genes_intro_bed_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/known_dom_genes_overlap_introgressed_regions.bed"

#-------------------
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
#genes_df <- read_excel(genes_fp)
genes_df <- read.csv(genes_fp)
# Load BED file with gene related info
bed_df <- read.delim(genes_intro_bed_fp, sep="\t", header=FALSE)
colnames(bed_df) <- c("Chrom", "Start", "End", "GenBankID")
# Split names, some are formatted as GenBankID_transcript_id
# Example: FJ974009.1-HORVU.MOREX.r3.7HG0637550.1 or GQ421469.1_cds_ACU68592.1_1
bed_df$GenBankID <- str_split_fixed(bed_df$GenBankID, '_', 2)[,1]
# Add gene symbols to bed df
intro.genes <- left_join(bed_df, genes_df, by=join_by("GenBankID"))
# Save to file
write_csv(intro.genes, file="known_dom_genes_introgressed_regions.csv", quote="none", col_names=TRUE)

# Write a BED file of the 50k SNPs in revised Table S2
genes_no_disease <- genes_df %>% dplyr::filter(Resistance_Gene != "yes" | is.na(Resistance_Gene))
#genes_no_disease %>% dplyr::select("SNP molecular marker", "Chrom.", "Position (bp)")
#bed_no_disease <- data.frame(Chrom=paste("chr", genes_no_disease$Chrom., sep=""), Start=as.numeric(genes_no_disease$`Position (bp)`)-1, End=as.numeric(genes_no_disease$`Position (bp)`), SNP=genes_no_disease$`SNP molecular marker`)
bed_no_disease <- data.frame(Chrom=genes_no_disease$chr_new, Start=as.numeric(genes_no_disease$pos_bp_new)-1, End=as.numeric(genes_no_disease$pos_bp_new), SNP=genes_no_disease$SNP.molecular.marker)
# Save to file the run bedtools intersect
write_delim(bed_no_disease, file="50k_markers_dom_genes-no_disease_resistance_genes.bed", delim="\t", na="", col_names=F)

# Use 50k_markers_dom_genes-no_disease_resistance_genes.bed file and bedtools to get 50k SNPs that overlap introgressed regions

# Load bed output
# BED file of 50k snps that overlap introgressed regions
snps_genes_intro_bed_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots_updated_chromPlots_our_pos/known_dom_genes_50k_snps_overlap_introgressed_regions.bed"
snps_bed_df <- read.delim(snps_genes_intro_bed_fp, sep="\t", header=FALSE)
colnames(snps_bed_df) <- c("Chrom", "Start", "End", "SNP_ID")
# Add gene symbols
intro.snps.genes <- left_join(snps_bed_df, genes_df, by=c("SNP_ID" = "SNP.molecular.marker")) %>%
    dplyr::select("Chrom.x", "Start.x", "End.x", "SNP_ID", "Locus.symbol")

#-------------------
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

#--------------------
# Define colors for groups
#color_breeding <- "#a70b0b"
color_breeding <- "#bf0e0e"
# color_wild <- "#bfdbf7"
color_wild <- "#64b5f6"
# color_hom_ref <- "#dee2e6"
"#faf9f9"
"#f8f9fa"
"#edede9"
"#dee2e6"
color_hom_ref <- "#edede9"
# color_genes <- "#353535"
#color_genes <- "#d6ccc2"
color_genes <- "#e9ecef"

#--------------------
# Prepare data for chromPlot
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
    dplyr::select(Chrom.x, Start.x, End.x, "Locus.symbol")
colnames(p.intro.genes) <- c("Chrom", "Start", "End", "ID")
# Add column of colors
p.intro.genes$Colors <- color_genes

# Prepare SNPs for plotting
colnames(intro.snps.genes) <- c("Chrom", "Start", "End", "SNP_ID", "ID")
# Make list of SNPs where we have gene intervals to exclude to avoid over plotting labels
p.intro.snps.genes <- intro.snps.genes[!(intro.snps.genes$SNP_ID %in% intro.genes$SNP.molecular.marker), ] %>%
    dplyr::select("Chrom", "Start", "End", "ID")
p.intro.snps.genes$Colors <- color_genes

p.intro.combo.genes <- rbind(p.intro.genes, p.intro.snps.genes)

setDT(p.intro.combo.genes)

for (samp in sample_names) {
    # Select sample
    curr_sample <- df_regions_gt[df_regions_gt$sample == samp, ] %>%
        dplyr::select(Chrom, Start, End, Name, Colors)
    # Pull genes that overlap introgressed segments only
    tmp_dom_curr_samp <- curr_sample[curr_sample$Name == "breeding", ]
    setDT(tmp_dom_curr_samp)
    setkey(tmp_dom_curr_samp, Chrom, Start, End)
    # Get region overlaps
    tmp_intro_genes <- foverlaps(p.intro.combo.genes, tmp_dom_curr_samp, type="any") %>%
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
                  stat=tmp_intro_genes, statCol="Value", statName="Value", noHist=TRUE,
                  cex=0.4)
    } else {
        chromPlot(gaps=chr_centromere, bands=curr_sample, figCols=7, title=samp)
    }
    #legend("bottomright", legend=c("breeding", "wild", "hom_reference"),
    #       col=c(color_breeding, color_wild, color_hom_ref), lty=1, lwd=5, cex=0.9)
    dev.off()
    # SVG
    # Generate and save plot
    out_svg <- paste0(samp, "_chrom_plot.svg")
    svg(out_svg, width = 7, height = 3)
    print(samp)
    if (nrow(tmp_intro_genes) > 1) {
        # There are genes that overlap introgressed regions for current sample
        chromPlot(gaps=chr_centromere, bands=curr_sample, figCols=7, title=samp,
                  stat=tmp_intro_genes, statCol="Value", statName="Value", noHist=TRUE,
                  cex=0.4)
    } else {
        chromPlot(gaps=chr_centromere, bands=curr_sample, figCols=7, title=samp)
    }
    dev.off()
}
