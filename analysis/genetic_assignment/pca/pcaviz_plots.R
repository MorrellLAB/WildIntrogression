library(PCAviz)
library(ggplot2)
library(magrittr)
library(cowplot)

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/plink_pca")

#--------------
# Exome capture PCA
excap.eigvec <- read.delim(file="dom_and_wbdc_exome_capture.pca.eigenvec", header = F, sep = " ")
excap.eigval <- read.delim(file="dom_and_wbdc_exome_capture.pca.eigenval", header = F, sep = " ")

plot(excap.eig)
