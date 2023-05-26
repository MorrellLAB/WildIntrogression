#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load java/openjdk-8_202
module load htslib/1.9
# Path to Beagle jar file
JAR_FILE="/panfs/jay/groups/9/morrellp/liux1299/Software/beagle.22Jul22.46e.jar"

# User provided input arguments
# Memory available
MEMORY="8"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed"

# Data parameters
GT_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz"
OUT_PREFIX="wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed"
# Plink MAP file with cM positions
#MAP_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.sorted.map"
MAP_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.excluded_problem_markers.map"
# Duplicate markers to exclude
EXCLUDE_MARKERS="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_duplicate_vars.dupvar"

# Phasing parameters
BURN_IN="3"
ITERATIONS="12"
PHASE_STATES="280"

# Imputation parameters
IMPUTE=true
IMPUTE_STATES="1600"
CLUSTER="0.005"
AP=true
GP=true

# General parameters
# Calculate effective population from nucleotide diversity estimates of theta
# Wild barley, theta estimate 0.008
Ne="400000"
EM=true
# Generate random seed and use here:
# echo $RANDOM
SEED="26015"

#---------------
mkdir -p ${OUT_DIR}

# Beagle phasing and imputation
java -Xmx${MEMORY}g -jar ${JAR_FILE} \
    gt=${GT_FILE} \
    out=${OUT_DIR}/${OUT_PREFIX} \
    map=${MAP_FILE} \
    excludemarkers=${EXCLUDE_MARKERS} \
    burnin=${BURN_IN} \
    iterations=${ITERATIONS} \
    phase-states=${PHASE_STATES} \
    impute=${IMPUTE} \
    imp-states=${IMPUTE_STATES} \
    cluster=${CLUSTER} \
    ap=${AP} \
    gp=${GP} \
    ne=${Ne} \
    em=${EM} \
    seed=${SEED}

tabix -p vcf --csi ${OUT_DIR}/${OUT_PREFIX}.vcf.gz