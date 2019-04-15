#!/bin/bash

set -e
set -o

#   This script performs Plink IBS/IBD analysis (IBS in this case because we don't have outgroup)

#   Usage: ./plink_ibs_chr.sh [VCF] [SNP_LIST] [MISSING] [FREQX_FILE] [OUT_DIR]

#   Where:
#   1) [VCF] is the full filepath to our VCF file
#   2) [SNP_LIST] is a list of full filepaths to files split into n SNP windows (i.e. 100 SNP windows)
#   3) [MISSING] is our missing genotype threshold (i.e. for 15% missing, use .15 here)
#   4) [FREQX_FILE] is a MAF file calculated by Plink (in another script) with .frqx extension
#   5) [OUT_DIR] is the full filepath to our output directory

#   Note: this script automatically generates output filenames based on prefix in SNP window files.
#   Output files will have the .genome file extension

#   Dependencies
module load plink/1.90b4.1

#   User provided arguments
VCF=$1
SNP_LIST=$2
MISSING=$3
FREQX_FILE=$4
OUT_DIR=$5

#   Generate from user provided arguments
SNP_LIST_ARR=($(cat ${SNP_LIST}))
#   Check if dir exists, if not make it
mkdir -p ${OUT_DIR}
mkdir -p ${OUT_DIR}/intermediates

#   Run Plink analysis
for s in ${SNP_LIST_ARR[@]}
do
    out_name=$(basename ${s} .txt)
    #   Now run IBS on each SNP window
    plink --vcf ${VCF} --extract ${s} --genome --geno ${MISSING} --read-freq ${FREQX_FILE} --allow-extra-chr --out ${OUT_DIR}/${out_name}
done

#   Do some cleanup and move .log and .nosex files to intermediates directory
mv ${OUT_DIR}/*.log ${OUT_DIR}/intermediates
mv ${OUT_DIR}/*.nosex ${OUT_DIR}/intermediates
