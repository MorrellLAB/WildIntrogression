#!/bin/bash

set -e
set -o pipefail

# Build the database: generate a transcript FASTA file.

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
GENE_PRED_FILE="$1"
REF_FASTA="$2"
OUT_DIR="$3"
OUT_PREFIX="$4"

#----------------
mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

# Run annovar script
retrieve_seq_from_fasta.pl --format refGene --seqfile ${REF_FASTA} ${GENE_PRED_FILE} --outfile ${OUT_PREFIX}_refGeneMrna.fa
