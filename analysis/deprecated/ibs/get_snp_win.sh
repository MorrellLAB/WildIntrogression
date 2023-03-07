#!/bin/bash

set -e
set -o pipefail

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
module load parallel
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/ibs

#   User provided arguments
SNP_LIST_BY_CHR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/snp_lists/all_chr_snp_file_list.txt"
WIN_SIZE="100"
SNP_SIZE="25"
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_ibs/snp_windows"

#-------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Generate SNP windows
parallel get_snp_win.py {} ${WIN_SIZE} ${SNP_SIZE} ${OUT_DIR} :::: ${SNP_LIST_BY_CHR}
