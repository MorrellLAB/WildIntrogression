#!/bin/bash

set -e
set -o pipefail

#   Dependencies
module load plink/1.90b6.10

#   User provided arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs"

#-----------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

#   Generated from user provided arguments
OUT_PREFIX=$(basename ${VCF} .vcf)

plink --vcf ${VCF} --allow-extra-chr --freqx --out ${OUT_DIR}/${OUT_PREFIX}
