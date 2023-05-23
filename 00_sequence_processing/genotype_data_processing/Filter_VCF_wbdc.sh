#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load gatk/4.1.2
module load htslib/1.9
# Functions used for counting number of sites in a VCF
source ~/GitHub/WildIntrogression/00_sequence_processing/count_num_sites_in_vcf.sh

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/wbdc_318_BOPA_morex_v3.vcf.gz"
OUT_PREFIX="wbdc_bopa_snps"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_introgression"

# Thresholds for filtering
# Proportion heterozygous genotypes threshold
HET_PROP="0.1"
# Max proportion missing
MAX_MISS="0.15"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${SCRATCH_DIR}

# Get the number of sites we are starting with
count_sites ${VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove sites that aren't polymorphic (minor allele count of 0) and unused alternate alleles.
echo "Removing sites that aren't polymorphic and unused alternate alleles..."
bcftools view --trim-alt-alleles --min-ac 1:nref ${VCF} -O z -o ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz
echo "Done removing sites that aren't polymorphic and unused alternate alleles."
# Get the number of sites left after filtering
count_sites ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with high missingness
echo "Removing sites that missingness exceeding threshold ${MAX_MISS}..."
bcftools filter -e "F_PASS(GT='mis') > ${MAX_MISS}" ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz -O z -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf.gz
MISS_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf.gz"
echo "Done removing sites with missingness exceeding threshold."
# Get the number of sites left after filtering and append to file
count_sites ${MISS_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter on heterozygosity (filter out highly heterozygous sites)
echo "Removing sites with more than ${HET_PROP} heterozygous genotypes..."
bcftools filter -i "COUNT(GT='het')/(N_SAMPLES-N_MISSING) < ${HET_PROP}" ${MISS_FILT_VCF} -O z -o ${OUT_DIR}/${OUT_PREFIX}.polymorphic.filt_miss_het.vcf.gz
echo "Done removing highly heterozygous sites."
HET_FILT_VCF="${OUT_DIR}/${OUT_PREFIX}.polymorphic.filt_miss_het.vcf.gz"
tabix -p vcf --csi ${HET_FILT_VCF}
# Get the number of sites left after filtering and append to file
count_sites ${HET_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
