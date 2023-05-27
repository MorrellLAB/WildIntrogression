#!/bin/bash

set -e
set -o pipefail

# Generate the genePred file from the GFF3 file

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
GFF_FILE="$1"
OUT_DIR="$2"
OUT_PREFIX="$3"

#-----------
mkdir -p ${OUT_DIR}

# Go into output dir
# Note: Get error when providing full output path directory instead of going into working dir
cd ${OUT_DIR}
# Generate the genePred file from the GFF3 file
gff3ToGenePred ${GFF_FILE} ${OUT_PREFIX}_refGene0.txt
# Sort genePred file
sort -V -k1,1 ${OUT_DIR}/${OUT_PREFIX}_refGene0.txt > ${OUT_DIR}/${OUT_PREFIX}_refGene0_sorted.txt

# Add a random first field to the genePred file using nl (this add the line number to each line)
nl ${OUT_DIR}/${OUT_PREFIX}_refGene0_sorted.txt > ${OUT_DIR}/${OUT_PREFIX}_refGene.txt

# Cleanup intermediate files
rm ${OUT_DIR}/${OUT_PREFIX}_refGene0.txt
rm ${OUT_DIR}/${OUT_PREFIX}_refGene0_sorted.txt
