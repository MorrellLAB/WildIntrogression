#!/bin/env bash

set -e
set -o pipefail

#   This script takes in Plink IBS/IBD output file (with .genome extension) and converts all spaces to tabs

# Dependencies
module load parallel/20210822

#   User provided arguments
# A list of the .genome files output from plink IBS/IBD analysis
GENOME_FILE_LIST="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/ibs_output/all_ibs_genome_list.txt"
# Full filepath to our output directory
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/ibs_output/reformatted_ibs_output"

#--------------------
mkdir -p ${OUT_DIR}

function reformat_ibs_file() {
    local genome_file="$1"
    local out_dir="$2"
    prefix=$(basename ${genome_file} .genome)
    #   Replace all spaces with tabs so output file is tab delimited
    awk -v OFS="\t" '$1=$1' ${genome_file} > ${out_dir}/${prefix}_reformat_genome.txt
}

export -f reformat_ibs_file

parallel --verbose reformat_ibs_file {} ${OUT_DIR} :::: ${GENOME_FILE_LIST}
