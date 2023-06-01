#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
module load htslib/1.9
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/genomics_general
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/genomics_general/VCF_processing

# User provided input arguments
VCF="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf.gz"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/abbababa"
OUT_PREFIX="dom_and_wild_morex_v3"
# Other settings
GENO_FORMAT="phased"
# Window size in bp
WIN_SIZE=200000000
MIN_SITES=100
STEP_SIZE=200000000
N_THREADS=6
MIN_PROP=0.5

#----------------------
# Convert VCF to the required .geno format
parseVCF.py -i ${VCF} | bgzip > ${OUT_DIR}/${OUT_PREFIX}.geno.gz

# Calculate ABBA BABA
ABBABABAwindows.py \
    -g ${OUT_DIR}/${OUT_PREFIX}.geno.gz \
    -o ${OUT_DIR}/${OUT_PREFIX}.csv \
    -f ${GENO_FORMAT} \
    -w ${WIN_SIZE} \
    -m ${MIN_SITES} \
    -s ${STEP_SIZE} \
    -P1 A \
    -P2 B \
    -P3 C \
    -O D \
    -T ${N_THREADS} \
    --minData ${MIN_PROP} \
    --popsFile pops.txt \
    --writeFailedWindows \
    --polarize &
