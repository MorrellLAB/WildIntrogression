#!/bin/bash

set -e
set -o pipefail

set -x # For debugging

# User provided input arguments
# Full path to a list of lists (to utilize GNU parallel and job arrays)
# List of lists here should only include FASTA files that had an alignment
#   (i.e., .fa and .tree files were generated) in the previous align step
FASTA_LIST_OF_LISTS="$1"
# Full path to the config file
CONFIG_FILE="$2"
# Full path to a list of MSA_Output directories that contain *.fa and *.tree files
MSA_DIR_LIST="$3"
# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script
SUBS_DIR="$4"
# List of subs names only that intersect with primary transcripts
# See script: intersect_primary_transcripts_and_subs.sh
#   Outputs the file: primary_transcript_intersect_names_only.txt
SUBS_NAMES_LIST="$5"
# Full path to output directory
OUT_DIR="$6"
# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT="$7"
# Full path to where we want to store the log files output from parallel
LOG_FILE_DIR="$8"

#------------------------------
function check_filepaths() {
    local curr_file="$1"
    # If file doesn't exist, return error and exit
    if ! [[ -f ${curr_file} ]]
    then
        echo "File doesn't exist: ${curr_file}, exiting..."
        exit 2
    fi
}

export -f check_filepaths

# Check if out dir exist, if not make them
mkdir -p ${OUT_DIR} ${LOG_FILE_DIR} ${OUT_DIR}/all_predict_log

# # Preparation steps
# # Generate list of substitutions files
# find ${SUBS_DIR} -name "*.subs" | sort -V > ${OUT_DIR}/subs_list.txt

# # Intersect primary transcripts and substitutions files
# if [ -f ${OUT_DIR}/primary_transcript_subs_list.txt ]
# then
#     echo "Primary transcript subs file list exists, proceeding with current list..."
# else
#     echo "Intersecting primary transcript names and substitutions files"
#     for pt in $(cat ${PRIMARY_TRANSCRIPTS})
#     do
#         if grep -wq ${pt} ${OUT_DIR}/subs_list.txt
#         then
#             # Primary transcript name is present in subs_list.txt
#             grep -w ${pt} ${OUT_DIR}/subs_list.txt >> ${OUT_DIR}/primary_transcript_subs_list.txt
#             # Also save just the primary transcript name due to naming scheme
#             #   (e.g., primary transcript is phvul.001g001100.1 and subs file is phvul.001g001100.1.v2.1.subs)
#             echo ${pt} >> ${OUT_DIR}/primary_transcript_intersect_names_only.txt
#         fi
#     done
# fi

# Build array of substutions files from primary transcripts only
#PT_SUBS_ARR=($(cat ${OUT_DIR}/primary_transcript_intersect_names_only.txt))
PT_SUBS_ARR=($(cat ${SUBS_NAMES_LIST}))
# Build array of MSA_output subdirectories containing .fa and .tree files
# Each task array index will be one current subdirectory (e.g., hvulgare_cds_list-000)
MSA_DIR_ARR=($(cat ${MSA_DIR_LIST}))

# Determine maximum array limit
MAX_ARRAY_LIMIT=$[${#MSA_DIR_ARR[@]} - 1]
echo "Maximum array limit is ${MAX_ARRAY_LIMIT}."

# Get current MSA_output subdirectory we are processing in current task array
CURR_MSA_DIR=$(echo ${MSA_DIR_ARR[${SLURM_ARRAY_TASK_ID}]})
echo "Current subdirectory we are processing in task array index ${SLURM_ARRAY_TASK_ID}: ${CURR_MSA_DIR}"
CURR_DIR_ID=$(basename ${CURR_MSA_DIR})
CURR_FASTA_LIST=$(grep -w ${CURR_DIR_ID} ${FASTA_LIST_OF_LISTS})

# Build array of MSA_output .fa and .tree files from current subdirectory (e.g., hvulgare_cds_list-000)
MSA_FASTA_ARR=()
for i in $(find ${CURR_MSA_DIR} -name "*.fasta" | sort -V)
do
    transcript_prefix=$(basename ${i} .fasta)
    # Match current fasta with primary transcript and subs files
    # Check if PT_SUBS_ARR contains current transcript fasta prefix, if so add to array for processing
    if [[ "${PT_SUBS_ARR[@]}" =~ "${transcript_prefix}" ]]
    then
        # Check to make sure .fasta and .tree files exist
        check_filepaths "${CURR_MSA_DIR}/${transcript_prefix}.fasta"
        # Add to array
        # Current fasta prefix is present in PT_SUBS_ARR
        MSA_FASTA_ARR+=("${CURR_MSA_DIR}/${transcript_prefix}.fasta")
    fi
done

function predict_sub() {
    local bad_mut_script="$1"
    local config_file="$2"
    local fasta_list="$3"
    local msa_fasta="$4"
    local curr_msa_dir="$5"
    local subs_list="$6"
    local out_dir="$7"
    # Get name of MSA_output fasta file
    msa_name=$(basename ${msa_fasta} .fasta)
    # Get msa_tree file
    msa_tree="${curr_msa_dir}/${msa_name}.tree"
    tree_name=$(basename ${msa_tree} .tree)

    # Check that msa_fasta and msa_tree names match
    if [[ ${msa_name} != ${tree_name} ]]
    then
        echo "There is a mismatch between the MSA fasta ${msa_name} and tree file ${tree_name}, exiting..."
        exit 1
    fi

    # Get current list prefix
    curr_list_prefix=$(basename ${fasta_list} .txt)
    # Get fasta file that corresponds with MSA fasta file
    fasta_file=$(grep -w ${msa_name} ${fasta_list})
    # Get subs filename
    subs_file=$(grep -w ${msa_name} ${subs_list})
    subdir_name="${curr_list_prefix}"
    # Check if out subdirectory exists, if not make it
    mkdir -p ${out_dir}/${subdir_name} ${out_dir}/all_predict_log/${subdir_name} ${out_dir}/${subdir_name}/problematic_predictions
    
    # Check if fasta file is in the list before if variable is empty
    if [ -z "$fasta_file" ]; then
        # Variable is empty
        # Look for fasta file in list before
        curr_list_num=$(echo $curr_list_prefix | cut -d'-' -f 2)
        if [ $curr_list_num != "000" ]; then
            # Look for fasta file in one list before
            # Remove leading zeroes
            list_num=$(echo $curr_list_num | sed 's/^0*//')
            # If current list number is less than or equal to 100
            if [ ${curr_list_num} == "100" ]; then
                # Add leading zero
                # Edge case when list_num=100 and previous list should be 099
                tmp_list_num=$(expr $list_num - 1)
                #prev_list_num=$(seq -f '%03g' ${tmp_list_num} ${tmp_list_num})
                prev_list_num=$(echo ${tmp_list_num} | sed 's/^/0/')
            else
                prev_list_num=$(expr $list_num - 1)
            fi
            # Use fasta list from one list before
            # Example:
            # Current list: /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/final_lists/hvulgare_cds_list-033.txt
            # One list before: /panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/final_lists/hvulgare_cds_list-032.txt
            prev_fasta_list=$(echo $fasta_list | sed -e "s/${list_num}.txt/${prev_list_num}.txt/")
            # Get fasta file that corresponds with MSA fasta file
            fasta_file=$(grep -w ${msa_name} ${prev_fasta_list})
        fi
    fi
    # Check that all files exist
    check_filepaths ${config_file}
    check_filepaths ${fasta_file}
    check_filepaths ${msa_fasta}
    check_filepaths ${msa_tree}
    check_filepaths ${subs_file}
    # echo "Config file: ${config_file}"
    # echo "Fasta file: ${fasta_file}"
    # echo "MSA fasta: ${msa_fasta}"
    # echo "MSA tree: ${msa_tree}"
    # echo "Subs file: ${subs_file}"

    # Predict substitutions
    # Redirect only stdout to log file and NOT stderr
    echo "Start predict step: ${out_dir}/${subdir_name}."
    set -x #Start debugging
    python ${bad_mut_script} predict \
        -c ${config_file} \
        -f ${fasta_file} \
        -a ${msa_fasta} \
        -r ${msa_tree} \
        -s ${subs_file} \
        -o ${out_dir}/${subdir_name} \
        1> ${out_dir}/all_predict_log/${subdir_name}/${msa_name}_predict.log
    set +x # End debugging
    # Stricter error checking so when predict file doesn't get written, the job doesn't
    #   return a zero exit status
    # Currently (2021-01-19), the predict output files should have the naming scheme: ${msa_name}_Predictions.txt
    #   Barley Example: HORVU0Hr1G000020.5_Predictions.txt
    if [[ -f ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt ]]
    then
        # Predict output file got written
        echo "Predict output file was written: ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt"
        echo "Checking predict output file..."
        # More in depth check of predict file in cases where file was written but error still occured.
        #   e.g., "Error: HyPhy killed by signal 15", "Function call stack", etc.
        if grep -wq "Error" ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt
        then
            # Move problematic prediction to subdirectory for easy troubleshooting
            mv ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt ${out_dir}/${subdir_name}/problematic_predictions
            echo "The following prediction resulted in an error: ${out_dir}/${subdir_name}/problematic_predictions/${msa_name}_Predictions.txt"
            echo "Please troubleshoot."
        fi
    else
        # Predict output file failed to get written to file, exiting
        echo "The following predict output file didn't successfully get written to a file: ${out_dir}/${subdir_name}/${msa_name}_Predictions.txt"
        echo "Exiting..."
        exit 2
    fi
}

export -f predict_sub

# Parallelize across array of paired MSA_output fasta and tree files
# Keep job log for parallel processes so upon resubmitting job, parallel can just re-run
#   samples that don't have an exit status of 0.
parallel --resume-failed --joblog ${LOG_FILE_DIR}/bad_mut_predict.sh.${SLURM_ARRAY_TASK_ID}.log predict_sub ${BAD_MUT_SCRIPT} ${CONFIG_FILE} ${CURR_FASTA_LIST} {} ${CURR_MSA_DIR} ${OUT_DIR}/primary_transcript_subs_list.txt ${OUT_DIR} ::: ${MSA_FASTA_ARR[@]}
