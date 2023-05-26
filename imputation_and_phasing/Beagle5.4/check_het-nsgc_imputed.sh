#!/bin/bash

# Dependencies
module load bcftools/1.10.2
module load texlive/20131202

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/vcf_summary"
# Arguments for bcftools stats
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"

SCRIPT_BCFTOOLS_STATS=~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering/run_bcftools_stats.sh

#---------------
source ${SCRIPT_BCFTOOLS_STATS}

# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR}

# Generate bcftools stats
generate_stats ${VCF} ${REF} ${OUT_DIR}
