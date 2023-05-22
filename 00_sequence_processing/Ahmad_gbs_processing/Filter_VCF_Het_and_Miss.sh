#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --tmp=4gb
#SBATCH -t 00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load htslib/1.9
# Functions used for counting number of sites in a VCF
source ~/GitHub/WildIntrogression/00_sequence_processing/count_num_sites_in_vcf.sh

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3/final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf.gz"
OUT_PREFIX="WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3"
# Thresholds for filtering
# Proportion heterozygous genotypes threshold
HET_PROP="0.1"
# Max proportion missing
MAX_MISS="0.20"
# Temporary directory (mainly to avoid going over storage quota)
SCRATCH_DIR="/scratch.global/liux1299/temp_introgression"

#----------------
mkdir -p ${OUT_DIR}

# Get the number of sites we are starting with
count_sites ${VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with high missingness
echo "Removing sites that missingness exceeding threshold ${MAX_MISS}..."
bcftools filter -e "F_PASS(GT='mis') > ${MAX_MISS}" ${VCF} -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf
MISS_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf"
echo "Done removing sites with missingness exceeding threshold."
# Get the number of sites left after filtering and append to file
count_sites ${MISS_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter on heterozygosity (filter out highly heterozygous sites)
echo "Removing sites with more than ${HET_PROP} heterozygous genotypes..."
bcftools filter -i "COUNT(GT='het')/(N_SAMPLES-N_MISSING) < ${HET_PROP}" ${MISS_FILT_VCF} -O z -o ${OUT_DIR}/${OUT_PREFIX}.filt_miss_het.vcf.gz
echo "Done removing highly heterozygous sites."
HET_FILT_VCF="${OUT_DIR}/${OUT_PREFIX}.filt_miss_het.vcf.gz"
tabix -p vcf --csi ${HET_FILT_VCF}
# Get the number of sites left after filtering and append to file
count_sites ${HET_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
