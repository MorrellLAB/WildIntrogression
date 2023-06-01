#!/usr/bin/env Rscript

# Import functions from another script
source("~/Software/genomics_general/jackknife.R")

allele_freq_fp <- "~/Dropbox/Projects/Wild_Introgression/Analyses/ABBA_BABA/nsgc_and_wbdc_geno/wild_and_dom.Hmurinum.derFreq.tsv.gz"
out_prefix <- "wild_and_dom.Hmurinum"

setwd("~/Dropbox/Projects/Wild_Introgression/Analyses/ABBA_BABA/nsgc_and_wbdc_geno")

# Case1: Ask whether there is evidence of introgression between P2 (wild_introgressed)
#   and P3 (domesticated). P1 (wild) will be the wild population
#   outgroup will be H. murinum (outgroup)
P1 <- "wild"
P2 <- "wild_introgressed"
P3 <- "domesticated"

#-----------------------------
# Genome-wide ABBA BABA
# Compute ABBA BABA proportions at each site, then
#   use these to compute the D statistic
D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA, na.rm=T) - sum(BABA, na.rm=T)) / (sum(ABBA, na.rm=T) + sum(BABA, na.rm=T))
}

# Read in data
freq_table <- read.table(allele_freq_fp, header=T, as.is=T)

nrow(freq_table)
head(freq_table)

# Compute D
D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])
# D varies from -1 to 1
# A positive D statistic would indicate an excess of ABBA over BABA
#   This suggests that P3 shares more genetic variation with P2 than with P1
print(paste("D =", round(D,4)))

# Use block-jackknife to test for a consistent genome-wide signal
block_indices <- get.block.indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)
print(paste("Genome divided into", n_blocks, "blocks."))

# Run the block jackknife procedure to compute the standard deviation of D
# D_sd <- get_jackknife_sd(block_indices=block_indices,
#                          FUN=D.stat,
#                          freq_table[,P1], freq_table[,P2], freq_table[,P3])
D_jackknife <- block.jackknife(block_indices=block_indices,
                               FUN=D.stat,
                               freq_table[,P1], freq_table[,P2], freq_table[,P3])

print(paste("D standard deviation = ", round(D_jackknife$mean,4)))

# From this unbiased estimate of the standard deviation of D, we can compute the standard error and
#   the Z score to test whether D deviates significantly from zero.
D_Z <- D_jackknife$mean / D_jackknife$standard_error
# D Z score 5.12, significant
print(paste("D Z score = ", round(D_Z,3)))

# Create table and save to file
dstat.out <- as.data.frame(D_jackknife)
dstat.out$D_zscore <- D_Z
write.table(dstat.out, file=paste0(out_prefix, ".D_and_Zscore.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#-----------------------------
# Estimate f, the admixture proportion
# Compare the excess of ABBA over BABA sites
# Function to compute f
# f.stat <- function(p1, p2, p3a, p3b) {
#   ABBA_numerator <- (1 - p1) * p2 * p3a
#   BABA_numerator <- p1 * (1 - p2) * p3a
#   
#   ABBA_denominator <- (1 - p1) * p3b * p3a
#   BABA_denominator <- p1 * (1 - p3b) * p3a
#   
#   (sum(ABBA_numerator) - sum(BABA_numerator)) /
#     (sum(ABBA_denominator) - sum(BABA_denominator))
# }
# 
# # Choose P3a and P3b and estimate f
# P3a <- "cyd_chi"
# P3b <- "cyd_zel"
# 
# f <- f.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])
# 
# print(paste("Admixture proportion = ", round(f,4)))
# 
# # Use block jackknife to estimate the standard deviation of f and obtain a confidence interval
# f_sd <- get_jackknife_sd(block_indices=block_indices,
#                          FUN=f.stat,
#                          freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])
# 
# # The 95% confidence interval is the mean +/- ~1.96 standard errors
# f_err <- f_sd/sqrt(n_blocks)
# 
# f_CI_lower <- f - 1.96*f_err
# f_CI_upper <- f + 1.96*f_err
# 
# print(paste("95% confidence interval of f =", round(f_CI_lower,4), round(f_CI_upper,4)))

#-----------------------------
# Chromosomal ABBA BABA
# Do all chromosomes show evidence of introgression?

# Identify all chromosome names present in the dataset
chrom_names <- unique(freq_table$scaffold)
chrom_indices <- lapply(chrom_names, function(chrom) which(freq_table$scaffold == chrom))
names(chrom_indices) <- chrom_names

# Check how many SNPs we have per chromosome
sapply(chrom_indices, length)

# Compute a D value for each chromosome
D_by_chrom <- sapply(chrom_names,
                     function(chrom) D.stat(freq_table[chrom_indices[[chrom]], P1],
                                            freq_table[chrom_indices[[chrom]], P2],
                                            freq_table[chrom_indices[[chrom]], P3]))

# Apply the jackknife to determine whether D differs significantly from zero for each chromosome
block_indices_by_chrom <- sapply(chrom_names,
                                 function(chrom) get.block.indices(block_size=1e6,
                                                                   positions=freq_table$position[freq_table$scaffold==chrom]),
                                 simplify=FALSE)

# Check the number of blocks per chromosome and the number of SNPs per block per chromosome
sapply(block_indices_by_chrom, length)
lapply(block_indices_by_chrom, sapply, length)

# Use jackknife to compute the Z scores for D for each chromosome
D_jackknife_by_chrom <- sapply(chrom_names,
                               function(chrom) block.jackknife(block_indices=block_indices_by_chrom[[chrom]],
                                                               FUN=D.stat,
                                                               freq_table[chrom_indices[[chrom]], P1],
                                                               freq_table[chrom_indices[[chrom]], P2],
                                                               freq_table[chrom_indices[[chrom]], P3]))

D_jackknife_by_chrom <- as.data.frame(t(D_jackknife_by_chrom))
D_jackknife_by_chrom$Z <- as.numeric(D_jackknife_by_chrom$mean) / as.numeric(D_jackknife_by_chrom$standard_error)
D_jackknife_by_chrom

# Create table and save to file
dstat.bychr <- data.frame(chr=rownames(D_jackknife_by_chrom),
                          mean=unlist(D_jackknife_by_chrom$mean),
                          variance=unlist(D_jackknife_by_chrom$variance),
                          standard_deviation=unlist(D_jackknife_by_chrom$standard_deviation),
                          standard_error=unlist(D_jackknife_by_chrom$standard_error),
                          Z=unlist(D_jackknife_by_chrom$Z))
write.table(dstat.bychr, file=paste0(out_prefix, ".D_and_Zscore.per_chr.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
