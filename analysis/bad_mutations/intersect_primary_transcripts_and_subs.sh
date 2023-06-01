#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
SUBS_DIR="$1"
OUT_DIR="$2"
PRIMARY_TRANSCRIPTS="$3"

# Check if out dir exist, if not make them
mkdir -p ${OUT_DIR}

# Preparation steps
# Generate list of substitutions files
find ${SUBS_DIR} -name "*.subs" | sort -V > ${OUT_DIR}/subs_list.txt

# Intersect primary transcripts and substitutions files
if [ -f ${OUT_DIR}/primary_transcript_subs_list.txt ]
then
    echo "Primary transcript subs file list exists, proceeding with current list..."
else
    echo "Intersecting primary transcript names and substitutions files"
    for pt in $(cat ${PRIMARY_TRANSCRIPTS})
    do
        if grep -wq ${pt} ${OUT_DIR}/subs_list.txt
        then
            # Primary transcript name is present in subs_list.txt
            grep -w ${pt} ${OUT_DIR}/subs_list.txt >> ${OUT_DIR}/primary_transcript_subs_list.txt
            # Also save just the primary transcript name due to naming scheme
            #   (e.g., primary transcript is phvul.001g001100.1 and subs file is phvul.001g001100.1.v2.1.subs)
            echo ${pt} >> ${OUT_DIR}/primary_transcript_intersect_names_only.txt
        fi
    done
fi
