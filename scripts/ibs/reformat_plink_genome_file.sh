#!/bin/env bash

#   This script takes in Plink IBS/IBD output file (with .genome extension) and converts all spaces to tabs

#   Usage: ./reformat_plink_genome_file.sh [GENOME_FILE] [OUT_DIR]

#   Where:
#   1) [GENOME_FILE] is the .genome file output form plink IBS/IBD analysis
#   2) [OUT_DIR] is the full filepath to our output directory

#   User provided arguments
GENOME_FILE=$1
OUT_DIR=$2

PREFIX=$(basename ${GENOME_FILE} .genome)

#   Replace all spaces with tabs so output file is tab delimited
awk -v OFS="\t" '$1=$1' ${GENOME_FILE} > ${OUT_DIR}/${PREFIX}_reformat_genome.txt
