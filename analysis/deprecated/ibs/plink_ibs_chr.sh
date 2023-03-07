#!/bin/bash

set -e
set -o pipefail

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
module load plink/1.90b6.10

#   User provided arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf"
SNP_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/snp_windows/all_chr_snp_windows_lists.txt"
MISSING="0.15"
FREQX_FILE="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/merged_domesticated_and_wbdc_318_morex_v3.frqx"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/ibs_output"

#--------------------------
# Make output directory
mkdir -p ${OUT_DIR} ${OUT_DIR}/intermediates

#   Generate from user provided arguments
SNP_LIST_ARR=($(cat ${SNP_LIST}))

#   Run Plink analysis
for s in ${SNP_LIST_ARR[@]}
do
    out_name=$(basename ${s} .txt)
    #   Now run IBS on each SNP window
    plink --vcf ${VCF} --extract ${s} --genome --geno ${MISSING} --read-freq ${FREQX_FILE} --allow-extra-chr --missing --out ${OUT_DIR}/${out_name}
done

#   Do some cleanup and move .log and .nosex files to intermediates directory
mv ${OUT_DIR}/*.log ${OUT_DIR}/intermediates
mv ${OUT_DIR}/*.nosex ${OUT_DIR}/intermediates
