#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=22gb
#SBATCH --tmp=6gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# Dependencies
module load admixture_ML/1.3.0
module load parallel/20210822

# User provided input arguments
PLINK_BED="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data_morex_v3/wbdc_bopa_snps_morex_v3_wPopInfo.pruned.bed"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/Admixture_wild_only_pruned"
# Min and Max K to try
MIN_K="2"
MAX_K="20"
# Number of replicate runs
NUM_REP_RUNS="10"

# Script that stores all the steps we need to run Admixture
ADMIXTURE_SCRIPT=~/GitHub/WildIntrogression/analysis/genetic_assignment/admixture.sh

#--------------
# Each replicate run with be submitted as a separate job through Slurm job arrays
# Within each job (one replicate), MIN_K to MAX_K values will be run in GNU parallel
# Prepare array of replicate runs
REP_RUNS_ARR=($(seq 1 ${NUM_REP_RUNS}))

# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#REP_RUNS_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current task array for current replicate run
CURR_REP_RUN=${REP_RUNS_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Current replicate run we are processing in task array index ${SLURM_ARRAY_TASK_ID}: run ${CURR_REP_RUN}"

${ADMIXTURE_SCRIPT} ${PLINK_BED} ${OUT_DIR}/run${CURR_REP_RUN} ${MIN_K} ${MAX_K}

# Run admixture with replicate runs in parallel
parallel --verbose ${ADMIXTURE_SCRIPT} ${PLINK_BED} ${OUT_DIR}/run{} ${MIN_K} ${MAX_K} ::: $(seq 1 ${NUM_REP_RUNS})
