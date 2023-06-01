#!/bin/bash

# Prepare .geno.gz file for ABBA BABA test
# See here for details: https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing
# This is a log of commands run

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
module load python2/2.7.12_anaconda4.2
GENOMICS_GENERAL="/panfs/roc/groups/9/morrellp/liux1299/Software/genomics_general/VCF_processing"

# User provided input arguments
VCF_SAMP="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf.gz"
VCF_OUTGROUP="/scratch.global/liux1299/barley_outgroups/Variant_Recalibrator/murinum_snps_final_pseudo.vcf.gz"
OUT_PREFIX="dom-wild-Hmurinum_morex_v3_pseudo"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/abbababa"

#-------------------
# Remove the following site from the REGIONS_FILE
# There is discordant REF/ALT alleles, so we'll exclude it
# Site: chr6H:42347342
OUT_NAME=$(basename ${VCF_SAMP} .vcf.gz)
VCF_DIR=$(dirname ${VCF_SAMP})
zgrep -vw "42347342" ${VCF_SAMP} > ${VCF_DIR}/${OUT_NAME}_concordant_sites.vcf
# Bgzip and index
bgzip ${VCF_DIR}/${OUT_NAME}_concordant_sites.vcf
tabix --csi -p vcf ${VCF_DIR}/${OUT_NAME}_concordant_sites.vcf.gz
# Set regions file
REGIONS_FILE="${VCF_DIR}/${OUT_NAME}_concordant_sites.vcf.gz"

# Merge vcf containing samples with vcf containing outgroup
bcftools merge -m all --regions-file ${REGIONS_FILE} -O v -o ${OUT_DIR}/${OUT_PREFIX}.vcf ${VCF_SAMP} ${VCF_OUTGROUP}

# Convert to ABBA BABA .geno.gz format
python ${GENOMICS_GENERAL}/parseVCF.py -i ${VCF_DIR}/${OUT_NAME}_concordant_sites.vcf.gz -i ~/scratch/barley_outgroups/Variant_Recalibrator/murinum_snps_final_pseudo.vcf.gz | bgzip > dom-wild-Hmurinum_morex_v3_pseudo.geno.gz
