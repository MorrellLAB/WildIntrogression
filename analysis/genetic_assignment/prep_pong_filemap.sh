#!/bin/bash

set -e
set -o pipefail



# Prepare filemap for pong
bed_prefix=$(basename ${PLINK_BED} .bed)
# Output to same directory as .Q files so it keeps things simple
# otherwise column3 would need to be path relative to the filemap
echo $(seq ${MIN_K} ${MAX_K}) | tr ' ' '\n' | awk -v file_prefix="${bed_prefix}" '{ print "k"$1 "\t" $1 "\t" file_prefix"."$1".Q" }' > ${OUT_DIR}/pong_filemap-${bed_prefix}.txt
