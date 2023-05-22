#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 60:00:00
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
module load python3/3.7.4_anaconda2019.10
module load vcftools_ML/0.1.16
module load bedtools/2.29.2
module unload R/3.4.4-tiff # One of the packages above also load an older version of R that will mess with downstream plotting

# User provided input arguments
VCF="/scratch.global/liux1299/temp_introgression/dom_and_wild_snps_polymorphic.vcf.gz"
OUT_PREFIX="dom_and_wild_snps"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_introgression"
# Thresholds for filtering
# Proportion heterozygous genotypes threshold
HET_PROP="0.1"
# Min and Max DP per sample threshold
MIN_DP="5"
MAX_DP="158"
# Max proportion missing per site
MAX_MISS="0.20"
# Quality cutoff
QUAL_CUTOFF="30"
# GQ cutoff per sample
GQ_CUTOFF="3"
# Allele balance filter, minimum and maximum cutoff for heterozygous genotypes
MIN_AB="0.40"
MAX_AB="0.60"
# Cutoff for filtering fake polymorphic sites
# These are sites that are actually monomorphic but show up as polymorphic
#   due to missingness and how programs handle that
# We'll filter out sites where the proportion is less than expected for a singleton
#   in the population. The minimum should be a singleton.
# To get this cutoff use: 1 / (number of samples * 2)
# Here, we have 614 total samples in the VCF: 1/1228 = 0.0008143322 (round down)
FAKE_POLY="0.0008"

# Uncallable BED file includes: REF has stretches of N's, repeat annotations, and high copy regions
UNCALLABLE_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/uncallable_regions/morex_v3_combined_uncallable.nochrUn.bed"

# Regions covered by exome capture (covered by >=50 reads)
# See Kono et al. 2019 for methods
# See Github repo for how this was generated:
# https://github.com/MorrellLAB/captured_50x_BED/tree/master/morex_v3
CAP50X_BED="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/captured_50x_morex_v3_partsRef_nochrUn.bed"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${SCRATCH_DIR} ${OUT_DIR}/vcf_summary

function check_filter_stringency() {
    local vcf_file="$1"
    local num_sites="$2"
    if [[ ${num_sites} == 0 ]]; then
        echo "No sites left after filtering. Try using less stringent criteria. File with no sites is: ${vcf_file}. Exiting..." >&2
        exit 8 # If not sites left, error out with message
    fi
}

export -f check_filter_stringency

function count_sites() {
    local vcf="$1"
    local log_file="$2"
    # Get the number of sites left after filtering
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        num_sites=$(zgrep -v "#" ${vcf} | wc -l)
    else
        # We are working with uncompressed vcf
        num_sites=$(grep -v "#" ${vcf} | wc -l)
    fi
    # Append the number of sites remaining to file
    printf "${vcf}\t${num_sites}\n" >> ${log_file}
    check_filter_stringency ${vcf} ${num_sites}
}

export -f count_sites

# Get the number of sites left for starting VCF and append to file
count_sites ${VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with GQ below cutoff, DP between minDP and maxDP, and by allelic balance
# Note: This script sets genotypes to missing if the GQ score is below the cutoff,
#   no sites actually get filtered out at this stage. Sites will get filtered out when we filter on missingness.
#   The same occurs for DP cutoffs, we filter out sites where depth is too low or too high
echo "Removing sites below GQ threshold ${GQ_CUTOFF}..."
echo "Also removing genotypes where per sample DP < ${MIN_DP} or per sample DP > ${MAX_DP}..."
echo "Filtering by min and max allelic balance where deviation for heterozygotes is ${MIN_AB} and ${MAX_AB}"
# Use bcftools to set genotypes to missing based on cutoffs
# Set GT to missing for het that fail allelic balance threshold
# AB defined the same way as Pedersen et al. 2021: alt/(ref+alt)
# Use bcftools to set genotypes to missing based on cutoffs
#   -i in +setGT means if GT ann meet condition, set to missing
#   FMT/AD[0:1] means first sample, second AD value
bcftools +setGT ${VCF} -- -t q -n "." -i "(GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])<${MIN_AB}) | (GT='het' & (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1])>${MAX_AB})" | bcftools +setGT - -- -t q -n "." -i "FMT/DP<${MIN_DP} | FMT/DP > ${MAX_DP} | FMT/GQ<${GQ_CUTOFF}" > ${SCRATCH_DIR}/${OUT_PREFIX}.AB_DP_GQ_filt.vcf
echo "Done setting sites below GQ threshold to missing."
echo "Done removing sites where per sample depth is either too low or too high."
AB_GQ_DP_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}.AB_DP_GQ_filt.vcf"

# Filter out sites with high missingness
echo "Removing sites that missingness exceeding threshold ${MAX_MISS}..."
bcftools filter -e "F_PASS(GT='mis') > ${MAX_MISS}" ${AB_GQ_DP_FILT_VCF} -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf
MISS_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtMISS.vcf"
echo "Done removing sites with missingness exceeding threshold."
# Get the number of sites left after filtering and append to file
count_sites ${MISS_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter on heterozygosity (filter out highly heterozygous sites)
echo "Removing sites with more than ${HET_PROP} heterozygous genotypes..."
bcftools filter -i "COUNT(GT='het')/(N_SAMPLES-N_MISSING) < ${HET_PROP}" ${MISS_FILT_VCF} -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtered_het.vcf
echo "Done removing highly heterozygous sites."
HET_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtered_het.vcf"
# Get the number of sites left after filtering and append to file
count_sites ${HET_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out sites with quality below cutoff
echo "Removing sites below quality threshold ${QUAL_CUTOFF}..."
bcftools filter -e "QUAL < ${QUAL_CUTOFF}" ${HET_FILT_VCF} -O z -o ${SCRATCH_DIR}/${OUT_PREFIX}_filtQUAL.vcf.gz
echo "Done removing sites below quality thresold."
QUAL_FILT_VCF="${SCRATCH_DIR}/${OUT_PREFIX}_filtQUAL.vcf.gz"
# Index vcf
tabix -p vcf ${QUAL_FILT_VCF}
# Get the number of sites left after filtering and append to file
count_sites ${QUAL_FILT_VCF} ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Filter out fake polymorphic sites
# These are sites that are actually monomorphic but got labeled as polymorphic
# This is just an extra check to be safe
echo "Removing fake polymorphic sites below MAF threshold ${FAKE_POLY}..."
#vcftools --gzvcf ${QUAL_FILT_VCF} --maf ${FAKE_POLY} --recode --recode-INFO-all --out ${OUT_DIR}/${OUT_PREFIX}_final
vcftools --gzvcf ${QUAL_FILT_VCF} --maf ${FAKE_POLY} --recode --recode-INFO-all --out ${OUT_DIR}/${OUT_PREFIX}_noFakePoly
echo "Done removing fake polymorphic sites."
# Rename output and compress
mv ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.recode.vcf ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf
bgzip ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf.gz

# Filter to only biallelic sites
echo "Filter to only biallelic sites..."
bcftools view -m2 -M2 ${OUT_DIR}/${OUT_PREFIX}_noFakePoly.vcf.gz -O z -o "${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz"
# Index VCF
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove variants that overlap with repeat annotated regions, high copy regions, and that overlap with stretches of Ns, and overlap with high diversity regions in morex-sample2
bedtools intersect -wa -v -header -a ${OUT_DIR}/${OUT_PREFIX}_biallelic.vcf.gz -b ${UNCALLABLE_BED} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Include only sites covered in exome capture
bedtools intersect -header -wa -a ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.vcf.gz -b ${CAP50X_BED} | bgzip > ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.cap50x.vcf.gz
tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.cap50x.vcf.gz
# Get the number of sites left after filtering and append to file
count_sites ${OUT_DIR}/${OUT_PREFIX}_biallelic.callable.cap50x.vcf.gz ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
