#!/usr/bin/env Rscript
#   Chaochih Liu - May 11, 2018

library(ggplot2)
library(scales)

readData <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE
    )
    return(df)
}

readPericentromere <- function(filename) {
    df <- read.delim(
        file = filename,
        header = FALSE,
        sep = " "
    )
    colnames(df) <- c("chr", "start", "end", "peri")
    return(df)
}

main <- function() {
    #   Take in command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    fp <- "/Users/chaochih/Projects/Wild_Introgression/Analyses/IBS/Results/pair_for_poster/CIho10420_WBDC172_introgression_sorted_merged_Int.txt"
    pericentromere.fp <- "/Users/chaochih/Dropbox/My_Posters/PEQG_2018/pericentromeres.txt"
    
    #   Read in file
    df <- readData(filename = fp)
    pericent.df <- readPericentromere(filename = pericentromere.fp)
    
    #   Key for pseudomolecular parts positions
    chr1H_part1 <- 312837513
    chr1H_part2 <- 245697919
    chr2H_part1 <- 393532674
    chr2H_part2 <- 374542350
    chr3H_part1 <- 394310633
    chr3H_part2 <- 305400481
    chr4H_part1 <- 355061206
    chr4H_part2 <- 291998952
    chr5H_part1 <- 380865482
    chr5H_part2 <- 289164678
    chr6H_part1 <- 294822070
    chr6H_part2 <- 288558443
    chr7H_part1 <- 325797516
    chr7H_part2 <- 331426484
    #   Whole chromosome size
    chr1.w <- chr1H_part1 + chr1H_part2
    chr2.w <- chr2H_part1 + chr2H_part2
    chr3.w <- chr3H_part1 + chr3H_part2
    chr4.w <- chr4H_part1 + chr4H_part2
    chr5.w <- chr5H_part1 + chr5H_part2
    chr6.w <- chr6H_part1 + chr6H_part2
    chr7.w <- chr7H_part1 + chr7H_part2
    #   Make chromosome sizes into data frame
    chrom_sizes <- data.frame(
        chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"),
        size = c(chr1.w, chr2.w, chr3.w, chr4.w, chr5.w, chr6.w, chr7.w)
    )
    
    chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7")
    chrom_key <- setNames(object = as.character(c(1:7)), nm = chrom_order)
    chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
    
    #   Convert chromosomes in each column to ordered factor
    chrom_sizes[["chr"]] <- factor(x = chrom_sizes[["chr"]], levels = chrom_order)
    df[["chr"]] <- factor(x = df[["chr"]], levels = chrom_order)
    pericent.df[["chr"]] <- factor(x = pericent.df[["chr"]], levels = chrom_order)
    
    #   Do the plotting
    ggplot(data = chrom_sizes) +
        #   base rectangles for chromosomes with numeric value for each chr on x-axis
        geom_rect(
            aes(
                xmin = as.numeric(chr) - 0.2,
                xmax = as.numeric(chr) + 0.2,
                ymax = size,
                ymin = 0
            ),
            colour = "gray80",
            fill = "gray95"
        ) +
        coord_flip() + # rotate plot 90 degrees
        #   black and white color theme
        theme(
            axis.text.x = element_text(colour = "black"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()
        ) +
        #   Give appearance of discrete axis with chr labels
        scale_x_discrete(labels = c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"), name = "Chromosome", limits = names(chrom_key)) +
        #   Add bands for pericentromeres
        geom_rect(
            data = pericent.df,
            aes(
                xmin = as.numeric(chr) - 0.20,
                xmax = as.numeric(chr) + 0.20,
                ymax = end,
                ymin = start
            ),
            colour = adjustcolor(col = "gray35", alpha.f = 0.5),
            fill = adjustcolor(col = "gray35", alpha.f = 0.5)
        ) +
        # scale_fill_manual(values) = group.colors +
        #   Add bands for introgressed regions
        geom_rect(
            data = df,
            aes(
                xmin = as.numeric(chr) - 0.20,
                xmax = as.numeric(chr) + 0.20,
                ymax = PhysPos_End,
                ymin = PhysPos_Start
            ),
            colour = adjustcolor(col = "skyblue", alpha.f = 0.6),
            fill = adjustcolor(col = "skyblue", alpha.f = 0.6)
        ) +
        ggtitle(paste0(unique(df$Ind1), " and ", unique(df$Ind2), " Introgressed Regions")) +
        #   Suppress Scientific notation on y-axis
        scale_y_continuous(labels = c(0, 200, 400, 600, 800)) +
        #scale_y_continuous(labels = comma) +
        ylab("Physical Position (Mb)")
}
