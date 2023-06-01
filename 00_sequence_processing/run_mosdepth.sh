#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=12gb
#SBATCH --tmp=10gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -e
set -o pipefail

# This script calculates per sample coverage utilizing Slurm job arrays

# Dependencies
# Mosdepth version 0.3.1
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/mosdepth

# User provided input arguments
BAM_LIST="/scratch.global/liux1299/temp_bam_introgression/all_bam_files_list.txt"
# Window size in bp for coverage estimation
WIN_SIZE="100"
REF_FASTA="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/coverage_exome"

#------------------
mkdir -p ${OUT_DIR}

# Prepare array for Slurm job array
BAM_ARR=($(cat ${BAM_LIST}))
# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#BAM_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get the current BAM file we are processing
CURR_BAM=${BAM_ARR[${SLURM_ARRAY_TASK_ID}]}
echo "Currently processing bam file: ${CURR_BAM}"

function run_mosdepth() {
    local bam_file="$1"
    local win_size="$2"
    local ref="$3"
    # Check if sample name contains substring
    if [[ "${bam_file}" == *"_phased_possorted_bam.bam"* ]]
    then
        # This is a 10x Genomics sample, generate clean sample name accordingly
        sample_prefix=$(basename ${bam_file} _phased_possorted_bam.bam)
    else
        sample_prefix=$(basename ${bam_file} .bam)
    fi
    set -x # For debugging
    mosdepth --by ${win_size} --threads 6 --fast-mode --fasta ${ref} --no-per-base ${sample_prefix} ${bam_file}
}

export -f run_mosdepth

# Go into out directory
cd ${OUT_DIR}
# Run mosdepth per sample
run_mosdepth ${CURR_BAM} ${WIN_SIZE} ${REF_FASTA}
