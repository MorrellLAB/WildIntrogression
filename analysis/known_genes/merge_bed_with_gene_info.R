#!/urs/bin/env Rscript

library(readxl)
library(dplyr)

setwd("~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes")

# Spreadsheet containing known domestication related genes
genes_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes/domestication-related genes12.14.2020Mike.xlsx"
# BED file of genes
genes_bed_fp <- "~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes/dom_genes_combined.gene_seq.cds.bed"

#--------------------
# Load gene related info
genes_df <- read_excel(genes_fp, sheet="main-updated-Mike2020", na=c("", "n/a")) %>%
  dplyr::select(GenBank_ID, Gene_Abbr, `Popular name`, `gene function`, `Gene role`, "Reference", `PubMed ID`, `Additional Notes`) %>%
  as.data.frame()
colnames(genes_df) <- c("GenBank_ID", "Gene_Abbr", "Popular_name", "gene_function", "Gene_role", "Reference", "PubMedID", "Additional_Notes")

# Load BED file
bed_df <- read.delim(genes_bed_fp, sep="\t", header=FALSE)
colnames(bed_df) <- c("Chrom", "Start", "End", "GenBank_ID")

# Split names, some are formatted as genbank_id-transcript_id
# Example: FJ974009.1-HORVU.MOREX.r3.7HG0637550.1
bed_df <- bed_df %>%
  separate(col=GenBank_ID, into=c("GenBank_ID", "Transcript_ID"), sep="-", fill="right")

# Add gene abbreviations to bed df
bed_genes <- full_join(bed_df, genes_df, by="GenBank_ID", multiple="all", relationship="many-to-many") %>%
  arrange(Chrom, Start, End, GenBank_ID, Gene_Abbr) %>%
  distinct()

# Get total count of unique gene names
length(unique(bed_genes$Gene_Abbr))

# Save to file
write_csv(x=bed_genes, file="domestication-related_genes.morex_v3.csv")
