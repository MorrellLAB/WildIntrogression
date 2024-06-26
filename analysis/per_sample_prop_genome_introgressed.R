#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(readxl)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/wild-introgressed-visible_pheno_and_Azerbaijan")
options(scipen=999)

fp <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v3/Table S5 - wbdc_likely_introgressed_segments-breeding.csv"
fp_passport <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v3/Table S1.xlsx"

# Morex v3 genome size in bp
morex_v3_size_bp <- 4225605719
# Mbp
morex_v3_size_mbp <- morex_v3_size_bp / 1000000

#-----------
df <- read.csv(fp)

df_passport <- read_excel(fp_passport) %>% dplyr::filter(grepl("WBDC", Accession_ID))

# Samples that show introgressed phenotype
vis_pheno <- df_passport %>% filter(Visible_Introgressed_Phenotype == "yes")
vis_pheno_samples <- vis_pheno$Accession_ID
# Pull wild-introgressed samples that show a visible phenotype
df_tracts_vis_pheno <- df[df$sample %in% vis_pheno_samples, ]
# Sum per sample introgressed segments
df_prop_vis_pheno <- df_tracts_vis_pheno %>%
    dplyr::select(sample, bp_length, Mbp_length) %>%
    group_by(sample) %>%
    dplyr::summarise(num_introgressed_segments=n(),
                     total_bp_length=sum(bp_length),
                     total_Mbp_length=sum(Mbp_length),
                     prop_genome_introgressed=sum(bp_length)/morex_v3_size_bp)
df_prop_vis_pheno$description <- "Visible Phenotype"

# Azerbaijan samples
azerbaijan_df <- df_passport %>% filter(Country_of_Origin == "Azerbaijan")
azerbaijan_samples <- azerbaijan_df$Accession_ID
# Pull Azerbaijan wild-introgressed samples
df_tracts_azerbaijan <- df[df$sample %in% azerbaijan_samples, ]
# Sum per sample introgressed segments
df_prop_azerbaijan <- df_tracts_azerbaijan %>%
    dplyr::select(sample, bp_length, Mbp_length) %>%
    group_by(sample) %>%
    dplyr::summarise(num_introgressed_segments=n(),
                     total_bp_length=sum(bp_length),
                     total_Mbp_length=sum(Mbp_length),
                     prop_genome_introgressed=sum(bp_length)/morex_v3_size_bp)
df_prop_azerbaijan$description <- "Azerbaijan"

# Wild-introgressed samples not with visible phenotype or from Azerbaijan
df_remaining <- df[!(df$sample %in% c(vis_pheno_samples, azerbaijan_samples)), ]
# Sum per sample introgressed segments
df_prop_remaining <- df_remaining %>%
    dplyr::select(sample, bp_length, Mbp_length) %>%
    group_by(sample) %>%
    dplyr::summarise(num_introgressed_segments=n(),
                     total_bp_length=sum(bp_length),
                     total_Mbp_length=sum(Mbp_length),
                     prop_genome_introgressed=sum(bp_length)/morex_v3_size_bp)
df_prop_remaining$description <- NA

df_out <- rbind(df_prop_vis_pheno, df_prop_azerbaijan)
df_all <- rbind(df_prop_vis_pheno, df_prop_azerbaijan, df_prop_remaining) %>%
    arrange(desc(prop_genome_introgressed))

sink(file="wild_introgressed_6number_summary-all.txt")
summary(df_all$prop_genome_introgressed)
sink()

# Save to file
write_csv(df_out, file="wild_introgressed_per_sample_prop-vis_pheno_and_Azerbaijan.csv", quote="none", col_names=TRUE)
write_csv(df_all, file="wild_introgressed_per_sample_prop-all.csv", quote="none", col_names=TRUE, na="")
