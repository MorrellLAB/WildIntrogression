#!/usr/bin/env Rscript
#   Original code by Li Lei March 15, 2018
#   Found at: https://github.com/MorrellLAB/Env_Assoc/blob/master/script/SNP_density/9k_exon_exon_capture_density_test1.R
#   Re-written and modified by Chaochih Liu - August 22, 2018

#   This script makes a fancy plot of BOPA and exon SNPs

library(ggplot2)
library(scales)

#   Define all functions needed first
ReadExData <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE,
        sep = "\t"
    )
    #   If chromsome column is named "Chromosome",
    #   change to "Chr", else do nothing
    if (names(df)[1] == "Chromosome") {
        names(df)[1] <- "Chr"
    }
    return(df)
}

ReadVcf <- function(filename) {
    df <- read.table(
        file = filename,
        header = FALSE,
        fill = TRUE,
        na.strings = "NA"
    )
    #   Extract only columns of interest
    vcfSubset <- data.frame(
        Chr = df$V1,
        PhysPos = df$V2,
        SNP_id = df$V3
    )
    return(vcfSubset)
}

RemoveChrUn <- function(dataFile) {
    #   This expects "chrUn" to represent unanchored chromosomes
    #   and Chromosome column to be "Chr" (this part is converted
    #   in ReadExData function)
    df <- dataFile[dataFile$Chr != "chrUn", ]
    return(df)
}

#   Do the work
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    excapFp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/SNP_density/ExomeCaptureTargets_per_Mb.txt"
    arraySnpFp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/SNP_density/GAPIT_bio6_GWAS_results_physPos.csv"
    exonSnpFp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/SNP_density/sorted_11_10380_Only_landrace_biallelic_NAM_final.pos"
    outDir <- "/Users/chaochih/Downloads"
    outName <- "test.pdf"
    
    #   This will get plotted as a line graph
    excapFp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/SNP_density/ExomeCaptureTargets_per_Mb.txt"
    #   These will be plotted as red triangles at the top
    #   (i.e. markers from our samples)
    arraySnpFp <- "/Users/chaochih/Dropbox/Projects/Wild_Introgression/vcf/merged_dom_and_wbdc_318_BOPA_sorted_noChrUn.vcf"
    
    #   Read in files
    excap <- ReadExData(filename = excapFp)
    arraySnp <- ReadVcf(filename = arraySnpFp)
    exonSnp <- ReadVcf(filename = exonSnpFp)
    #exonSnp <- read.delim(file = exonSnpFp)
    
    #   Drop any unmapped chromosome
    excapNoUn <- RemoveChrUn(dataFile = excap)
    arraySnpNoUn <- RemoveChrUn(dataFile = arraySnp)
    arraySnpNoUn <- arraySnpNoUn[arraySnpNoUn$Chr != "1/1", ]
    exonSnpNoUn <- RemoveChrUn(dataFile = exonSnp)
    exonSnpNoUn <- exonSnpNoUn[exonSnpNoUn$Chr != "specificone", ]
    
    #   Create plot
    pdf(file = paste0(outDir, "/", outName), width = 10, height = 6)
    
    #   Plotted as vertical light blue lines
    ggplot(excapNoUn) +
        #geom_vline(aes(xintercept = as.numeric(PhysPos)/1000000), size = 0.02, alpha = 0.1, color = "#a6cee3") +
        #   Plotted as blue line graph
        geom_line(aes(x = (Start+End)/2000000, y = NExCap/10), size = 0.75, alpha = 0.7, color = "#1f78b4") +
        #   Plotted as red triangles
        geom_point(aes(x = as.numeric(PhysPos)/1000000, y = 15), data = arraySnpNoUn, shape = 17, size = 0.95, alpha = 0.45, color = "#ff4000") +
        facet_grid(Chr~.) +
        scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
        scale_x_continuous(limits = c(0, 725), breaks = seq(0, 725, by = 50)) +
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text.y = element_text(size = 10, colour = "black", angle = 0),
            axis.ticks.y = element_blank()) +
        labs(y = "", x = "Physical Position (Mb)")
    
    dev.off()
}
       
main() # Run the program