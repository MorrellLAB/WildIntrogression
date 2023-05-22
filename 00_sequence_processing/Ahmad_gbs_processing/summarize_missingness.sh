#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load vcftools_ML/0.1.16
module load R/4.0.4

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3/final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf.gz"
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3/vcf_summary"

# We are calling on some helper scripts from sequence_handling
#   and will need to define a few variables
seqhand="/panfs/jay/groups/9/morrellp/liux1299/sequence_handling"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${OUT_DIR}/Intermediates

function calc_missing() {
    local vcf="$1"
    local out_dir="$2"
    local out_prefix="$3"
    if [[ ${vcf} == *".gz"* ]]; then
    # Use gzip flag
    vcftools --gzvcf ${vcf} \
        --missing-indv \
        --out ${out_dir}/${out_prefix}_missingness
    # Missing data per site
    vcftools --gzvcf ${vcf} \
        --missing-site \
        --out ${out_dir}/${out_prefix}_missingness
    else
        # Uncompressed VCF
        vcftools --vcf ${vcf} \
            --missing-indv \
            --out ${out_dir}/${out_prefix}_missingness
        # Missing data per site
        vcftools --vcf ${vcf} \
            --missing-site \
            --out ${out_dir}/${out_prefix}_missingness
    fi
}

export -f calc_missing

VCF_PREFIX=$(basename ${VCF} .vcf.gz)

# VCF
calc_missing ${VCF} ${OUT_DIR} ${VCF_PREFIX}
# Visualize missingness
"${seqhand}/HelperScripts/graph_missingness.R" \
    ${OUT_DIR}/${VCF_PREFIX}_missingness.imiss \
    ${OUT_DIR}/${VCF_PREFIX}_missingness.lmiss \
    ${OUT_DIR}
