#!/bin/bash

set -e
set -o

#   Dependencies
module load plink/1.90b4.1

#   User provided arguments
VCF=$1
MISSING=$2
OUT_DIR=$3

#   Generated from user provided arguments
OUT_PREFIX=$(basename ${VCF} .vcf)

plink --vcf ${VCF} --allow-extra-chr --freqx --out ${OUT_DIR}/${OUT_PREFIX}
