#!/bin/bash

# Dependencies
module load vcflib_ML/1.0.0_rc2
module load python3/3.8.3_anaconda2020.07_mamba
module load java/openjdk-8_202
# Path to flare.jar file
JAR_FILE="/panfs/jay/groups/9/morrellp/liux1299/Software/flare.jar"

# User provided input arguments
# VCF file with phased reference genotypes
REF_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_merged.phased.imputed.no_missing.vcf.gz"
# VCF file with phased genotypes to be analyzed
GT_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/wild_likely_introgressed.polymorphic.filt_miss_het.phased.imputed.vcf.gz"
# PLINK map file with cM units
MAP="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_and_dom_snps.polymorphic.filt_miss_het.excluded_problem_markers.map"
# file with reference sample to panel map
REF_PANEL="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/FLARE/dom_and_wild_not_introgressed.ref.panel"
# Output file prefix
OUT_PREFIX="wild_likely_introgressed_geno"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/local_ancestry-flare"
# Duplicate markers to exclude
EXCLUDE_MARKERS="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_and_dom_duplicate_vars.dupvar"

# General parameters
# Is input data from a SNP array
ARRAY=true

# Generate random seed and use here:
# echo $RANDOM
SEED="28610"

#------------------
mkdir -p ${OUT_DIR}

java -jar ${JAR_FILE} \
    ref=${REF_VCF} \
    gt=${GT_VCF} \
    map=${MAP} \
    ref-panel=${REF_PANEL} \
    array=${ARRAY} \
    seed=${SEED} \
    excludemarkers=${EXCLUDE_MARKERS} \
    out=${OUT_DIR}/${OUT_PREFIX}.flare.out

# Convert output VCF to BED format using vcflib utility
zcat ${OUT_DIR}/${OUT_PREFIX}.flare.out.anc.vcf.gz | vcf2bed.py > ${OUT_DIR}/${OUT_PREFIX}.flare.out.anc.bed
