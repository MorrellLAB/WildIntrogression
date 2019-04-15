#!/usr/bin/env Rscript

#   Chaochih Liu - May 7, 2018

#   This script takes in a .genome file outputted from Plink IBS/IBD analysis and generates a histogram
#   and quantile of Z2 values. The purpose of this script is purely for data exploration to determine
#   a Z2 threshold to use, so the plots are not publication ready as is.

library(data.table)

readFile <- function(filename) {
    df <- fread(
        input = filename,
        header = TRUE
    )
    return(df)
}

main <- function() {
    #   Take in command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    #   For testing
    fp <- args[1]
    out.dir <- args[2]
    
    #   Read in data
    df <- readFile(filename = fp)
    
    bn <- basename(fp)
    out.prefix <- sub(
        pattern = ".genome",
        x = bn,
        replacement = "",
        ignore.case = TRUE
    )
    
    pdf(file = paste0(out.dir, "/", out.prefix, "_Z2_hist.pdf"))
    hist(
        x = df$Z2,
        xlim = c(0, 1),
        breaks = 100,
        xlab = "Plink IBS/IBD Z2 Value",
        main = paste0("Histogram of", out.prefix)
    )
    dev.off()
    
    q <- quantile(x = df$Z2, probs = seq(0, 1, 0.05))
    
    write.table(
        x = q,
        file = paste0(out.dir, "/", out.prefix, "_Z2_quantile.txt"),
        quote = FALSE,
        sep = "\t"
    )
}

main()
