#!/bin/bash

set -e
set -o pipefail

# This script compiles all predictions into a single file for easy downstream processing

# User provided input arguments
# Full path to a list of predict output directories that contain *_Predictions.txt files
#   This can be generated as follows from within the predict script ${OUT_DIR}:
#   cd /path/to/predict_output_${SAMPLE_NAME}
#   Generate the list of output directories
#   find $(pwd -P) -maxdepth 1 -type d -name "*list*" > predict_out_dir_list.txt
#   Get path to list of output directories
#   find $(pwd -P) -name predict_out_dir_list.txt
PREDICT_DIR_LIST="$1"
# Full path to the list of long substitutions .txt file
#   i.e., The long_substitutions.txt file output from the script VeP_to_Subs.py
LONG_SUBS_FILE="$2"
# Project name will be used as prefix for final combined report
PROJECT="$3"
# Where do we want to store our final compiled predictions report?
OUT_DIR="$4"
# Full path to the BAD_Mutations.py script
BAD_MUT_SCRIPT="$5"

#------------------------------
# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Prepare array
PREDICT_DIR_ARR=($(cat ${PREDICT_DIR_LIST}))

# Function compiles predictions for one subdirectory
function compile_predictions() {
    local bad_mut_script="$1"
    local curr_predict_subdir="$2"
    local long_subs_file="$3"
    echo "Compiling predictions for subdirectory: ${curr_predict_subdir}"
    python ${bad_mut_script} compile \
        -P ${curr_predict_subdir} \
        -S ${long_subs_file}
}

export -f compile_predictions

# Parallelize across predict subdirectories
parallel --verbose compile_predictions ${BAD_MUT_SCRIPT} {} ${LONG_SUBS_FILE} ::: ${PREDICT_DIR_ARR[@]}

# BAD_Mutations.py compile outputs combined reports in each subdirectory
# Now, compile reports in subdirectories into a single file
# First, get the base directory of the predict output subdirectories
#   Assume base directory is the same for all subdirectories
PREDICT_BASE_DIR=$(dirname ${PREDICT_DIR_ARR[0]})

# Make a list of Combined_Report.txt files in subdirectories to compile into single file
find ${PREDICT_BASE_DIR} -name "Combined_Report.txt" | sort -V > ${OUT_DIR}/${PROJECT}_predict_subdir_reports_list.txt
# Store in array
REPORTS_ARR=($(cat ${OUT_DIR}/${PROJECT}_predict_subdir_reports_list.txt))

# Prepare header line for final compiled reports file
head -n 1 ${REPORTS_ARR[0]} > ${OUT_DIR}/${PROJECT}_Combined_Report.txt

# Append combined reports excluding header line to final compiled reports file
for i in ${REPORTS_ARR[@]}
do
    grep -v "VariantID" ${i} >> ${OUT_DIR}/${PROJECT}_Combined_Report.txt
done
