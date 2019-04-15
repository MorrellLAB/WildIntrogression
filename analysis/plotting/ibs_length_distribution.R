#!/usr/bin/env Rscript

library(ggplot2)
library(beeswarm)
library(ggbeeswarm)
library(scales)

ReadData <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE,
        sep = "\t"
    )
    return(df)
}

main <- function() {
    #   Take in command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    dataFp <- args[1]
    outDir <- args[2]
    outFileName <- args[3] # include file extension
    
    #   Read in data
    ibsData <- ReadData(filename = dataFp)
    
    pdf(file = paste0(outDir, "/", outFileName), width = 8.40, height = 5.40)
    
    #   Make plot
    ggplot(ibsData, aes(x = factor(Chr), y = Int_Phys_Size/1000000)) +
        geom_violin(aes(fill = Chr)) +
        scale_fill_brewer() +
        scale_y_continuous(
            name = "IBS Interval Size (Mbp)",
            labels = c("0", "100", "200", "300", "400", "500", "600"),
            breaks = c(0, 100, 200, 300, 400, 500, 600),
            limits = c(0, 600)
        ) +
        scale_x_discrete(
            name = "Chromosome",
            labels = c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
        ) +
        labs(title = "Distribution of introgression intervals across chromosomes") +
        theme(plot.title = element_text(size = 16))
    
    dev.off()       
}

main()
