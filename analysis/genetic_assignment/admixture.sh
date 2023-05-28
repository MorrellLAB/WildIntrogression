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
    random_seed=${RANDOM}
    admixture --seed=${random_seed} ${plink_bed} ${kval} --haploid="*"
    admixture --seed=${random_seed} --cv ${plink_bed} ${kval} | tee log${kval}.out
}

export -f run_admixture

parallel --verbose run_admixture ${PLINK_BED} {} ::: $(seq ${MIN_K} ${MAX_K})
