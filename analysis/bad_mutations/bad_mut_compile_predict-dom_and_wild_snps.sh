#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=6gb
#SBATCH --tmp=2gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o bad_mut_compile_predictions.sh.%j.out
#SBATCH -e bad_mut_compile_predictions.sh.%j.err

set -e
set -o pipefail

# This script compiles all predictions in parallel into a single file for easy downstream processing.
#   BAD_Mutations.py will output the compiled report in the predict output subdirectories. These will then be compiled into a single report in ${OUT_DIR}.

# Usage: sbatch bad_mut_compile_predict-${PROJECT}.job

# Dependencies
module load parallel
module load python3/3.8.3_anaconda2020.07_mamba

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User provided input arguments
# Full path to a list of predict output directories that contain *_Predictions.txt files
#   This can be generated as follows from within the predict script ${OUT_DIR}:
#   cd /path/to/predict_output_${SAMPLE_NAME}
#   Generate the list of output directories
#       find $(pwd -P) -maxdepth 1 -type d -name "*list*" | sort -V > predict_out_dir_list.txt
#   Get path to list of output directories
#       realpath predict_out_dir_list.txt
PREDICT_DIR_LIST="/scratch.global/liux1299/bad_mutations/predict_output_dom_and_wild_snps/predict_out_dir_list.txt"

# Full path to the list of long substitutions .txt file
#   i.e., The long_substitutions.txt file either generate from output from the script VeP_to_Subs.py or Annovar
LONG_SUBS_FILE="/scratch.global/liux1299/bad_mutations/predict_output_dom_and_wild_snps/subs_list.txt"
# Project name will be used as prefix for final combined report
PROJECT="dom_and_wild_snps"

# Where do we want to store our final compiled predictions report?
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations/predictions"

# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/Software/BAD_Mutations/BAD_Mutations.py"

# Full path to the script bad_mut_compile_predictions.sh
COMPILE_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/bad_mutations/bad_mut_compile_predictions.sh"

#------------------------------
# Run compile predictions script
${COMPILE_SCRIPT} \
    ${PREDICT_DIR_LIST} \
    ${LONG_SUBS_FILE} \
    ${PROJECT} \
    ${OUT_DIR} \
    ${BAD_MUT_SCRIPT}
