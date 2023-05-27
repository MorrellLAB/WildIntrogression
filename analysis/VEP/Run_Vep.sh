#!/bin/bash

set -e
set -o pipefail

# Define function to run VeP, this will be called in other job scripts
function run_vep() {
    local vcf=$1
    local gff_file=$2
    local ref=$3
    local species=$4
    # Prepare out prefix from filename
    out_prefix=$(basename ${vcf} .vcf.gz)
    # Run VeP
    # Don't use --total_length flag, turning on messes up BAD_Mutations Vep_to_Subs.py script
    vep -i ${vcf} \
        --format vcf \
        --check_svs \
        --custom ${gff_file},,gff \
        --fasta ${ref} \
        --gff ${gff_file} \
        --output_file ${out_prefix} \
        --species ${species} \
        --verbose \
        --force \
        --warning_file ${out_prefix}.log
    # Rename txt file output
    mv ${out_prefix} ${out_prefix}.txt
}

export -f run_vep
