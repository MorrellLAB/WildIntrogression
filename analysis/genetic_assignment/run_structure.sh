#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# written by Li Lei, 20190416, adapted from Peter L. Morrel's script: https://github.com/pmorrell/Utilities/blob/master/Structure_replicates.sh
# Modified by Chaochih Liu to utilize Slurm job arrays and GNU parallel
#   to speed up processing multiple K values and reps

# Usage: sbatch --array=0-8 run_structure.sh
# Where: 0-8 corresponds to the number of array indices from KMIN to KMAX
#   Example: KMIN=2 and KMAX=10 generates an array where the maximum
#   array limit is 8 (array containing a sequence of values from 2 to 10).

set -e
set -o pipefail

# Dependencies
module load structure/2.3.3
module load parallel/20210822

# User provided input arguments
# Structure data file
STRCT_IN="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data/wbdc_318_BOPA_revised.strct_in.txt"
# Structure main parameters file
MAIN_PARAMS="/panfs/roc/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/genetic_assignment/318_wild_mainparams"
# Structure extra parameters
EXTRA_PARAMS="/panfs/roc/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/genetic_assignment/318_wild_extraparams"
# Output directory
OUT_DIR="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/Structure_wild_only"
OUT_PREFIX="wild_318"
REPS="3"
KMIN="2"
KMAX="10"

#------------------------
mkdir -p ${OUT_DIR}

# Prepare array of K values
K_ARR=($(seq ${KMIN} ${KMAX}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#K_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."
# Get K for current Slurm job array index
CURR_K=${K_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently running K of: ${CURR_K}"

# Prepare array of the number of replicates
REP_ARR=($(seq 1 ${REPS}))

# For current K, pick a random number of seconds from 1-10 to delay parallel jobs
# This is to make sure each K and each replicate has a unique seed
num_sec=$(seq 1 10 | shuf | head -n 1)
num_sleep_sec=$(seq 1 900 | shuf | head -n 1)
# Sleep for s seconds before starting parallel processing in case mutliple
#   Slurm job arrays start at the same time, this is to make sure we are using a unique seed
sleep ${num_sleep_sec}s
# Run current K for n reps in parallel
# Delay starting each job by n seconds so seed is unique for each replicate run
parallel --verbose --delay ${num_sec} structure -i ${STRCT_IN} -m ${MAIN_PARAMS} -e ${EXTRA_PARAMS} -K ${CURR_K} -o ${OUT_DIR}/${OUT_PREFIX}_K${CURR_K}_rep{} ::: ${REP_ARR[@]}
