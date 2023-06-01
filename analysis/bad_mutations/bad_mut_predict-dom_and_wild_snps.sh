#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10gb
#SBATCH --tmp=2gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o bad_mut_predict.sh.%A_%a.out
#SBATCH -e bad_mut_predict.sh.%A_%a.err

set -e
set -o pipefail

# This script stores all filepaths and calls on the script bad_mut_predict.sh

# Dependencies
module load parallel
module load python3/3.8.3_anaconda2020.07_mamba

# Activate conda environment
source activate /home/morrellp/liux1299/.conda/envs/bad_mutations

# User provided input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
# List of lists here should only include FASTA files that had an alignment
#	(i.e., .fa and .tree files were generated) in the previous align step
FASTA_LIST_OF_LISTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/final_lists/all_cds_hvulgare_list_of_lists.txt"

# Full path to the config file
CONFIG_FILE="/panfs/jay/groups/9/morrellp/gfrascar/bad_mutations_scripts/config.txt"

# Full path to a list of MSA_Output directories that contain *.fa and *.tree files
MSA_DIR_LIST="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/MSA_output/MSA_out_dir_list_of_dir.txt"

# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script or ANNOVAR_to_subs.py script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations/annovar_to_subs/dom_and_wild_snps"

# List of subs names only that intersect with primary transcripts
# See script: intersect_primary_transcripts_and_subs.sh
#   Outputs the file: primary_transcript_intersect_names_only.txt
SUBS_NAMES_LIST="/scratch.global/liux1299/bad_mutations/predict_output_dom_and_wild_snps/primary_transcript_intersect_names_only.txt"

# Sample name will be used as a prefix for outputs
SAMPLE_NAME="dom_and_wild_snps"

# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"

# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/Software/BAD_Mutations/BAD_Mutations.py"

# Full path to where we want to store the log files output from parallel
LOG_FILE_DIR="${OUT_DIR}/all_parallel_log_files"

# Script that stores predict_sub function
PREDICT_SCRIPT="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/bad_mutations/bad_mut_predict.sh"

#------------------------------
# Run predict script
${PREDICT_SCRIPT} ${FASTA_LIST_OF_LISTS} \
	${CONFIG_FILE} \
	${MSA_DIR_LIST} \
	${SUBS_DIR} \
	${SUBS_NAMES_LIST} \
	${OUT_DIR} \
	${BAD_MUT_SCRIPT} \
	${LOG_FILE_DIR}
