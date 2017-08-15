#!/usr/bin/env Rscript
# Connor Depies Aug, 15, 2017


# Takes 3 arguments and produces the mean of each column of a lower triangular matrix from a symmetric square matrix.
#   To Run: ./lowerTri_colMinMax.R [plink.dist] [outDir] [outFilename]

#   Required arguments:
#       1) plinkDistFile A file containing a square symmetric distance matrix.
#       2) outDir: where do we want our output files to go?
#           IMPORTANT Caveat: Make sure to remove the last "/" at the end of the filepath
#       3) outFilename: output file name excluding file extension
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
calcMean <- function(tri.matrix) {
    cat("Calculating mean of each column...", sep = "\n")
    #   Convert to data frame
    tri.matrix.df <- as.data.frame(tri.matrix)
    #   Calculate mean
    tmp.mean <- apply(
        X = tri.matrix.df,
        MARGIN = 2, # iterate over columns
        FUN = mean, # average
        na.rm = TRUE # remove NA's
    )
    #   Store mean values in data frame
    m.df <- data.frame(MEAN = tmp.mean)
    cat("Done with calculations and saving to data frame...", sep = "\n")
    return(m.df)
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
#   Stores arguments into a vector
plinkDistFile <- args[1] # 1) Plink.dist/plink.mibs to be turned into triangular matrix
outDirectory <- args[2] # 2) Where should we output our files?
outFilename <- args[3] # 3) Output filename (EXCLUDE file extension)
#   Read in file
plink.m <- readFile(filename = plinkDistFile)
# Extracts a lower triangular matrix from the square matrix
m.lt <- replaceUpperTri(square.matrix = plink.m) 
# Calculates mean
m.df <- calcMean(tri.matrix = m.lt)
# Outputs to file
writeToOut(out.df = m.df, outDir = outDirectory, outputName = outFilename) # Writes to output




