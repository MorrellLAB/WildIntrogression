#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(writexl)

old_fp <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v1/Table S1.xlsx"
full_passport_fp <- "~/Dropbox/Projects/Wild_Introgression/Steffenson_Comments/ALL_WBDC_PASSPORT_INFO_MAY-2023.xlsx"

# Introgressed sample lists
# Wild-introgressed identified in Fang et al 2014
fang_2014_fp <- "~/GitHub/WildIntrogression/analysis/trees/Fang_et_al_2014_wild-introgressed_list.txt"
this_study_fp <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v3/Table S3 - per_sample_counts_of_introgressed_regions.csv"

#-------------
df_old <- read_excel(old_fp, sheet=1)
df_passport <- read_excel(full_passport_fp, sheet=1)
df_passport <- df_passport[1:318, ]

df_fang_2014 <- read.delim(fang_2014_fp, header = F)
colnames(df_fang_2014)[1] <- "Accession_ID"
df_fang_2014$Identified_in_Previous_Study <- "Fang et al 2014"

df_this_study <- read.csv(this_study_fp) %>% select(sample)
colnames(df_this_study)[1] <- "Accession_ID"
df_this_study$Identified_in_This_Study <- "x"
df_this_study$Approaches <- "FLARE"

df_reordered <- df_old[, c("Accession_ID", "Country_of_Origin", "Longitude", "Latitude", "Accession_Type")]

df_dom <- df_reordered[!grepl("WBDC", df_reordered$Accession_ID), ]
df_old_wild <- df_reordered[grepl("WBDC", df_reordered$Accession_ID), ]

temp_sample_names <- df_passport[1:318, c(2, 4, 5, 17, 20, 27)]
# Check if sample names are all the same
unique(temp_sample_names[1] == temp_sample_names[2])
unique(temp_sample_names[1] == temp_sample_names[3])
unique(temp_sample_names[1] == temp_sample_names[4])
unique(temp_sample_names[1] == temp_sample_names[5])
unique(temp_sample_names[1] == temp_sample_names[6])
# All columns ahve the same sample names, pick one for merge

df_sub <- data.frame(Accession_ID=gsub("-", "", df_passport$`Original WBDC Accession to ICARDA Ig Number Assignment`), Country_of_Origin=df_passport$`Country of originc`, Longitude=df_passport$`Longitude lon_dd`, Latitude=df_passport$`Latitude lat_dd`)
df_sub$Accession_Type <- "wild"

df_updated <- rbind(df_sub, df_dom)

# Add wild-introgressed columns
temp_df <- left_join(df_updated, df_fang_2014, by="Accession_ID")
df_updated_introg <- left_join(temp_df, df_this_study, by="Accession_ID")

writexl::write_xlsx(df_updated_introg, path="~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v3/Table S1.xlsx")
