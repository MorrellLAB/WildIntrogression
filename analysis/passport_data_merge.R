#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)

old_fp <- "~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v1/Table S1.xlsx"
full_passport_fp <- "~/Dropbox/Projects/Wild_Introgression/Steffenson_Comments/ALL_WBDC_PASSPORT_INFO_MAY-2023.xlsx"

#-------------
df_old <- read_excel(old_fp, sheet=1)
df_passport <- read_excel(full_passport_fp, sheet=1)
df_passport <- df_passport[1:318, ]

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

write_csv(df_updated, file="~/Dropbox/Projects/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v3/Table S1.csv", na="")
