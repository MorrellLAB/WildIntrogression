#!/bin/bash

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
GENOMICS_GENERAL="/panfs/jay/groups/9/morrellp/liux1299/Software/genomics_general"

# User provided input arguments
GENO_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/abbababa/nsgc_and_wbdc_geno/dom-wild-Hmurinum_morex_v3.geno.gz"
POPS_FILE="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/abbababa/nsgc_and_wbdc_geno/dom-wild-Hmurinum.pops.txt"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/dxy"
OUT_PREFIX="dom-wild-Hmurinum_morex_v3"

WIN_SIZE="20000"
MIN_SITES_PER_WIN="50"

THREADS="5"
MIN_SITES_PER_WIN="200"

#--------
mkdir -p ${OUT_DIR}

python ${GENOMICS_GENERAL}/popgenWindows.py -w ${WIN_SIZE} -m ${MIN_SITES_PER_WIN} -g ${GENO_FILE} --outFile ${OUT_DIR}/${OUT_PREFIX}.win${WIN_SIZE}.minsites${MIN_SITES_PER_WIN}.csv.gz -f phased -T ${THREADS} -p wild -p wild_introgressed -p domesticated --popsFile ${POPS_FILE}
