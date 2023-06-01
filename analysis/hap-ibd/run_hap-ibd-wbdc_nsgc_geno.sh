#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=6gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load java/openjdk-8_202
# Path to flare.jar file
JAR_FILE="/panfs/jay/groups/9/morrellp/liux1299/Software/hap-ibd/test/hap-ibd.jar"

# User provided input arguments
# VCF file with phased genotypes to be analyzed
GT_VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz"
# PLINK map file with cM units
MAP="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_and_dom_snps.polymorphic.filt_miss_het.excluded_problem_markers.map"
# Output file prefix
OUT_PREFIX="dom_and_wild_with_introgressed_geno"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/hap_ibd"

# Max gap
# Try max intermarker distance from VCF
MAX_GAP="1000"
# Minimum number of markers for an IBS segment
MIN_MARKERS="20"

#------------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/minmarkers${MIN_MARKERS}_maxgap${MAX_GAP}

java -jar ${JAR_FILE} \
    gt=${GT_VCF} \
    map=${MAP} \
    out=${OUT_DIR}/minmarkers${MIN_MARKERS}_maxgap${MAX_GAP}/${OUT_PREFIX}.hap-ibd.out \
    max-gap=${MAX_GAP} \
    min-markers=${MIN_MARKERS}
