#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load parallel

#   User provided arguments
VCF="$1"
OUT_DIR="$2"

function split_by_chr() {
    local vcf="$1"
    local chrom="$2"
    local out_dir="$3"
    prefix=$(basename ${vcf} .vcf.gz)
    bcftools filter ${vcf} -r ${chrom} -O z -o ${out_dir}/${prefix}-${chrom}.vcf.gz
}

export -f split_by_chr

# Prepare array
chrom_arr=($(zgrep -v "#" ${VCF} | cut -f 1 | sort -uV))

# Split VCF by chromosome
parallel split_by_chr ${VCF} {} ${OUT_DIR} ::: ${chrom_arr[@]}
