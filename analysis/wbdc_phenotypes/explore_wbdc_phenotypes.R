#!/usr/bin/env Rscript

library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# User provided input arguments
fp_aug2007 <- "~/Dropbox/Projects/Wild_Introgression/Data/Phenotype_Data-Brian_Steffenson/Morrell_WBDC_morph_data_Aug26_07.xls"
fp_wbdc_2021 <- "~/Dropbox/Projects/Wild_Introgression/Data/Phenotype_Data-Brian_Steffenson/Unpublished_2021/THE_WBDC_FINAL.xlsx"
# Working directory - output plots and tables here
setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/wbdc_phenotypes")

#------------------
# Read in data and convert to long format
df_aug2007 <- read_xls(fp_aug2007, sheet="Data_for_analysis", na="-") %>%
  mutate_at(c('Plant_height', 'spike_length', 'flag_leaf_length', 'flag_leaf_width', 'awn_length', 'tiller_number', 'flag-1_leaf_length', 'flag-1_leaf_width', 'heading_days', 'Kernel_number'), as.numeric) %>%
  as.data.frame()
# Rename columns
colnames(df_aug2007) <- c('ID', "Habitat", 'Plant_height', 'Spike_length', 'Flag_leaf_length', 'Flag_leaf_width', 'Awn_length', 'Tiller_number', 'Flag1_leaf_length', 'Flag1_leaf_width', 'Heading_days', 'Kernel_number')

# Explore outliers
df_aug2007 %>% dplyr::select(ID, Flag_leaf_length) %>% arrange(desc(Flag_leaf_length)) %>% dplyr::top_n(5)
df_aug2007 %>% dplyr::select(ID, Flag_leaf_width) %>% arrange(desc(Flag_leaf_width)) %>% dplyr::top_n(5)
df_aug2007 %>% dplyr::select(ID, Tiller_number) %>% arrange(Tiller_number) %>% dplyr::top_n(-10)

# Subset trait measurements
# Flag leaf length
flag_leaf_length_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Flag_leaf_length, Flag1_leaf_length) %>%
  pivot_longer(cols=c('Flag_leaf_length', 'Flag1_leaf_length'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Flag leaf width
flag_leaf_width_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Flag_leaf_width, Flag1_leaf_width) %>%
  pivot_longer(cols=c('Flag_leaf_width', 'Flag1_leaf_width'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Spike length and awn length
spike_awn_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Spike_length, Awn_length) %>%
  pivot_longer(cols=c('Spike_length', 'Awn_length'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Plant height
plant_height_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Plant_height) %>%
  pivot_longer(cols=c('Plant_height'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Tiller number, Kernel number
tiller_kernel_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Tiller_number, Kernel_number) %>%
  pivot_longer(cols=c('Tiller_number', 'Kernel_number'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Heading days
heading_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Heading_days) %>%
  pivot_longer(cols=c('Heading_days'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()

# Leaf length plot
ggplot(flag_leaf_length_traits_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_beeswarm(color="#6b9080", cex=1.3, alpha=0.8) +
  theme_classic() +
  ylab("Length (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(flag_leaf_length_traits_07$Measurement))
# Save plot
ggsave("wbdc_aug07-flag_leaf_length.jpg", width=6, height=5, units="in", dpi=300)

# Leaf width plot
ggplot(flag_leaf_width_traits_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_quasirandom(color="#6b9080", alpha=0.8) +
  theme_classic() +
  ylab("Width (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(flag_leaf_width_traits_07$Measurement))
ggsave("wbdc_aug07-flag_leaf_width.jpg", width=6, height=5, units="in", dpi=300)

# Spike length and Awn length plot
ggplot(spike_awn_traits_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_quasirandom(color="#6b9080", alpha=0.8) +
  theme_classic() +
  ylab("Length (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(spike_awn_traits_07$Measurement))
ggsave("wbdc_aug07-spike_and_awn_length.jpg", width=6, height=5, units="in", dpi=300)

# Plant height plot
ggplot(plant_height_trait_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_quasirandom(color="#6b9080", alpha=0.8) +
  theme_classic() +
  ylab("Height (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(plant_height_trait_07$Measurement))
ggsave("wbdc_aug07-plant_height.jpg", width=6, height=5, units="in", dpi=300)

# Tiller and kernel number
ggplot(tiller_kernel_trait_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_quasirandom(color="#6b9080", alpha=0.8) +
  theme_classic() +
  ylab("Count") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(tiller_kernel_trait_07$Measurement))
ggsave("wbdc_aug07-tiller_and_kernel_number.jpg", width=6, height=5, units="in", dpi=300)

# Heading days
ggplot(heading_trait_07, aes(x=Trait, y=Measurement)) +
  geom_boxplot(outlier.shape=NA, color="#736f72") +
  geom_quasirandom(color="#6b9080", alpha=0.8) +
  theme_classic() +
  ylab("Days") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(heading_trait_07$Measurement))
ggsave("wbdc_aug07-heading_days.jpg", width=6, height=5, units="in", dpi=300)

#------------------
# Read in data
df_leaf1and2 <- readxlsx(fp_wbdc_2021, sheet="1st Leaf and 2nd Leaf")
