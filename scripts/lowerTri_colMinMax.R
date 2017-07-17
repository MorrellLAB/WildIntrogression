#!/usr/bin/env Rscript

#   This is a script that extracts a lower triangular matrix
#       from a square matrix (excludes diagonal) and calculates
#       the minimum and maximum of each column
#   Written by Chaochih Liu
#   July 17, 2017

#   To Run: ./lowerTri_colMinMax.R [plink.dist] [outDir] [outFilename]

#   Required arguments:
#       1) plink.dist: square matrix output from plink
#       2) outDir: where do we want our output files to go?
#           IMPORTANT Caveat: Make sure to remove the last "/" at the end of the filepath
#       3) outFilename: output file name excluding file extension

#   Function to read in plink square matrix
readFile <- function(filename) {
    #   read in file as data frame
    #   convert to matrix
    m <- as.matrix(read.table(
        file = filename,
        fill = TRUE
    ))
    cat("Done reading in file...", sep = "\n")
    return(m)
}

#   Function to extract lower triangular matrix for calculations
replaceUpperTri <- function(square.matrix) {
    cat("Extracting lower triangular matrix...", sep = "\n")
    #   We want to work with the lower matrix
    #   Replace upper triangular matrix (including diagonal) with "NA" values
    square.matrix[upper.tri(square.matrix, diag = TRUE)] <- NA
    #   Remove last column of all NA's before calculation
    #   This column should always be empty given a symmetric square matrix
    s.matrix <- square.matrix[, -ncol(square.matrix)]
    cat("Done extracting lower triangular matrix...", sep = "\n")
    return(s.matrix)
}

#   Function to calculate minimum and maximum of each column
#   Store in data frame
calcMinMax <- function(tri.matrix) {
    cat("Calculating minimum and maximum of each column...", sep = "\n")
    #   Convert to data frame
    tri.matrix.df <- as.data.frame(tri.matrix)
    #   Calculate min first
    tmp.min <- apply(
        X = tri.matrix.df,
        MARGIN = 2, # iterate over columns
        FUN = min, # minimum
        na.rm = TRUE # remove NA's
    )
    #   Then calculate max
    tmp.max <- apply(
        X = tri.matrix.df,
        MARGIN = 2, # iterate over columns
        FUN = max, # maximum
        na.rm = TRUE # remove NA's
    )
    #   Store min and max values in data frame
    mm.df <- data.frame(Minimum = tmp.min, Maximum = tmp.max)
    cat("Done with calculations and saving to data frame...", sep = "\n")
    return(mm.df)
}

#   Function to save data to output file
writeToOut <- function(out.df, outDir, outputName) {
    cat("Saving data to output files...", sep = "\n")
    #   Add file extension to output filename
    name <- paste0(outDir, "/", outputName, ".txt")
    write.table(
        x = out.df,
        file = name, # what is the full filepath and output filename?
        quote = FALSE,
        sep = "\t", # output file will be tab delimited
        eol = "\n",
        col.names = TRUE, # include column names
        row.names = FALSE
    )
    cat("Done.", sep = "\n")
}

#   Do the work
#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)
#   User provided arguments
plinkDistFile <- args[1] # 1) plink.dist square matrix
outDirectory <- args[2] # 2) Where should we output our files?
outFilename <- args[3] # 3) Output filename (EXCLUDE file extension)

#   Read in file
plink.m <- readFile(filename = plinkDistFile)
#   Extract lower triangular matrix
m.lt <- replaceUpperTri(square.matrix = plink.m)
#   Calculate minimum and maximum
mm.df <- calcMinMax(tri.matrix = m.lt)
#   Save data to output file
writeToOut(out.df = mm.df, outDir = outDirectory, outputName = outFilename)
