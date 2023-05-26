#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load admixture_ML/1.3.0

# User provided input arguments
PLINK_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data_morex_v3/wbdc_bopa_snps_morex_v3_wPopInfo.pruned.bed"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/Admixture_wild_only_pruned"
# Min and Max K to try
MIN_K="2"
MAX_K="20"

#--------------
mkdir -p ${OUT_DIR}
# Go into output directory, Admixture outputs files to current working dir
cd ${OUT_DIR}

# Run admixture for various K values
for kval in $(seq ${MIN_K} ${MAX_K})
do
    admixture ${PLINK_BED} ${kval} --haploid="*"
done

# Calculate the standard error so we can pick a K value
for kval in $(seq ${MIN_K} ${MAX_K})
do
    admixture --cv ${PLINK_BED} ${kval} | tee log${kval}.out
done
