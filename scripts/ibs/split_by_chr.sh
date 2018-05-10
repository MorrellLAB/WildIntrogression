#!/bin/bash

set -e
set -o

#   User provided arguments
VCF=$1
OUT_DIR=$2

function splitByChr() {
    local vcf=$1
    local out_dir=$2
    #   Skip vcf metadata ('##') and header ('#') lines
    #   Extract only CHROM and ID columns from VCF file
    #   Remove unknown chromosomes
    #   OFS is for output field separater
    grep -v "#" ${vcf} | awk -v OFS='\t' '{print $1,$3}' | grep -v "chrUn" > ${out_dir}/tmp_chr_and_snps_only.txt

    #   Split file by chromosome (7 chromosomes for barley)
    for i in $(seq 1 7)
    do
        grep "chr${i}H" ${out_dir}/tmp_chr_and_snps_only.txt | cut -f 2 > ${out_dir}/chr${i}_snpsList.txt
    done

    #   Cleanup tmp files
    rm ${out_dir}/tmp_chr_and_snps_only.txt
}

export -f splitByChr

#   Run functions and do the work
#   Extract SNP windows
splitByChr ${VCF} ${OUT_DIR}
