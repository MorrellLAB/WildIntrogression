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

# Identified introgressed samples
intro.samp <- c("WBDC 016", "WBDC 020", "WBDC 027", "WBDC 042", "WBDC 045", "WBDC 053", "WBDC 092", "WBDC 125", "WBDC 128", "WBDC 172", "WBDC 182", "WBDC 190", "WBDC 227", "WBDC 275", "WBDC 280", "WBDC 283", "WBDC 311", "WBDC 344")

#------------------
# Read in data and convert to long format
df_aug2007 <- read_xls(fp_aug2007, sheet="Data_for_analysis", na="-") %>%
  mutate_at(c('Plant_height', 'spike_length', 'flag_leaf_length', 'flag_leaf_width', 'awn_length', 'tiller_number', 'flag-1_leaf_length', 'flag-1_leaf_width', 'heading_days', 'Kernel_number'), as.numeric) %>%
  as.data.frame()
# Rename columns
colnames(df_aug2007) <- c('ID', "Habitat", 'Plant_height', 'Spike_length', 'Flag_leaf_length', 'Flag_leaf_width', 'Awn_length', 'Tiller_number', 'Flag1_leaf_length', 'Flag1_leaf_width', 'Heading_days', 'Kernel_number')

# Add introgressed sample labels
df_aug2007 <- df_aug2007 %>%
  mutate(
    Pop = case_when(
      ID %in% intro.samp ~ "introgressed",
      .default = "wild"
    )
  )

# Explore outliers
df_aug2007 %>% dplyr::select(ID, Pop, Flag_leaf_length) %>% arrange(desc(Flag_leaf_length)) %>% dplyr::top_n(5)
df_aug2007 %>% dplyr::select(ID, Pop, Flag_leaf_width) %>% arrange(desc(Flag_leaf_width)) %>% dplyr::top_n(5)
df_aug2007 %>% dplyr::select(ID, Pop, Tiller_number) %>% arrange(Tiller_number) %>% dplyr::top_n(-10)

# Subset trait measurements
# Flag leaf length
flag_leaf_length_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Flag_leaf_length, Flag1_leaf_length) %>%
  pivot_longer(cols=c('Flag_leaf_length', 'Flag1_leaf_length'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Flag leaf width
flag_leaf_width_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Flag_leaf_width, Flag1_leaf_width) %>%
  pivot_longer(cols=c('Flag_leaf_width', 'Flag1_leaf_width'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Spike length and awn length
spike_awn_traits_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Spike_length, Awn_length) %>%
  pivot_longer(cols=c('Spike_length', 'Awn_length'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Plant height
plant_height_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Plant_height) %>%
  pivot_longer(cols=c('Plant_height'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Tiller number, Kernel number
tiller_kernel_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Tiller_number, Kernel_number) %>%
  pivot_longer(cols=c('Tiller_number', 'Kernel_number'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()
# Heading days
heading_trait_07 <- df_aug2007 %>%
  dplyr::select(ID, Pop, Heading_days) %>%
  pivot_longer(cols=c('Heading_days'),
               names_to="Trait",
               values_to="Measurement") %>%
  as.data.frame()

# Leaf length plot
ggplot(flag_leaf_length_traits_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Length (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(flag_leaf_length_traits_07$Measurement))
# Save plot
ggsave("wbdc_aug07-flag_leaf_length.jpg", width=6, height=5, units="in", dpi=300)

# Leaf width plot
ggplot(flag_leaf_width_traits_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Width (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(flag_leaf_width_traits_07$Measurement))
ggsave("wbdc_aug07-flag_leaf_width.jpg", width=6, height=5, units="in", dpi=300)

# Spike length and Awn length plot
ggplot(spike_awn_traits_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Length (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(spike_awn_traits_07$Measurement))
ggsave("wbdc_aug07-spike_and_awn_length.jpg", width=6, height=5, units="in", dpi=300)

# Plant height plot
ggplot(plant_height_trait_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Height (mm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(plant_height_trait_07$Measurement))
ggsave("wbdc_aug07-plant_height.jpg", width=6, height=5, units="in", dpi=300)

# Tiller and kernel number
ggplot(tiller_kernel_trait_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Count") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(tiller_kernel_trait_07$Measurement))
ggsave("wbdc_aug07-tiller_and_kernel_number.jpg", width=6, height=5, units="in", dpi=300)

# Heading days
ggplot(heading_trait_07, aes(x=Trait, y=Measurement, fill=Pop)) +
  geom_violin(width=0.8, position=position_dodge(0.8), alpha=0.5) +
  geom_boxplot(width=0.1, alpha=0.8, position=position_dodge(0.8)) +
  scale_fill_manual(values=c("#cce3de", "#3c6e71")) +
  theme_classic() +
  ylab("Days") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  ylim(0, max(heading_trait_07$Measurement))
ggsave("wbdc_aug07-heading_days.jpg", width=6, height=5, units="in", dpi=300)

#------------------
# Statistical tests
# Two sample non-parametric t-test
flwidth.wild <- flag_leaf_width_traits_07[flag_leaf_width_traits_07$Pop == "wild", "Measurement"] %>% na.omit()
flwidth.introgressed <- flag_leaf_width_traits_07[flag_leaf_width_traits_07$Pop == "introgressed", "Measurement"] %>% na.omit()
shapiro.test(flwidth.wild)
shapiro.test(flwidth.introgressed)

flag_leaf_width_traits_07 %>%
  group_by(Pop) %>%
  dplyr::summarise(
    count=n(), median=median(Measurement, na.rm=TRUE), IQR=IQR(Measurement, na.rm=TRUE)
  )

wilcox.test(flwidth.introgressed, flwidth.wild)

var.f.flwidth <- var(flwidth.introgressed) / var(flwidth.wild)
pf(q=var.f.flwidth, df1=length(flwidth.introgressed)-1, df2=length(flwidth.wild)-1, lower.tail=F)

# Awn lengths
awn.wild <- spike_awn_traits_07[spike_awn_traits_07$Pop == "wild", "Measurement"] %>% na.omit()
awn.introgressed <- spike_awn_traits_07[spike_awn_traits_07$Pop == "introgressed", "Measurement"] %>% na.omit()

shapiro.test(awn.wild)
shapiro.test(awn.introgressed)

spike_awn_traits_07 %>%
  group_by(Pop) %>%
  dplyr::summarise(
    count=n(), median=median(Measurement, na.rm=TRUE), IQR=IQR(Measurement, na.rm=TRUE)
  )

#------------------
# Read in data
df_leaf1and2 <- readxlsx(fp_wbdc_2021, sheet="1st Leaf and 2nd Leaf")
