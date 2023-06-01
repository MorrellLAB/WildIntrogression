#!/bin/bash

# User provided input arguments
PLINK_BED="$1"
OUT_DIR="$2"
# Min and Max K to try
MIN_K="$3"
MAX_K="$4"

#--------------
mkdir -p ${OUT_DIR}
# Go into output directory, Admixture outputs files to current working dir
cd ${OUT_DIR}

# Run admixture for various K values
# Parallelize across K values
function run_admixture() {
    local plink_bed=$1
    local kval=$2
    local prefix=$3
    random_seed=${RANDOM}
    admixture --seed=${random_seed} ${plink_bed} ${kval} --haploid="*"
    admixture --cv ${plink_bed} ${kval} | tee log.${prefix}.${kval}.out
}

export -f run_admixture

prefix=$(basename ${PLINK_BED} .bed)

parallel --verbose run_admixture ${PLINK_BED} {} ${prefix} ::: $(seq ${MIN_K} ${MAX_K})

# Compile cross-validation error results across values of K
echo "# CV results" > ${OUT_DIR}/${prefix}.CV.txt
for ((K=${MIN_K};K<=${MAX_K};K++)); do
    awk -v K=$K '$1=="CV"{ print K,$4 }' ${OUT_DIR}/log.${prefix}.${K}.out >> ${OUT_DIR}/${prefix}.CV.txt
done
