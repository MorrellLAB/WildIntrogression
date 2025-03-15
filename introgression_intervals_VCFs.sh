#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Original script Peter L. Morrell - St. Paul, MN 
# This script is designed to extract the genotyping data for a set of samples from a VCF file
# The script will generate a VCF file for each sample in the list

set -e
set -o pipefail

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

# Load necessary modules
# module load bcftools

# Define the input files
BED=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/introgression_intervals.bed
# Genotyping data file
VCF=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/dom_and_wild_compare.vcf.gz
# All samples to be included except the query (hybrid) sample
SAMPLE_LIST=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/BOPA_cultivated_samples.txt

OUTPUT_DIR="/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/Results"

log "   -> Generating list of samples for comparisons"
SAMPLES=($(cut -f4 "${BED}" | sort -u))

log "   -> Removing single SNP intervals"
BED_INTERVALS=$(awk -F'\t' '$5 != 1' "${BED}")

create_sample_vcf() {
    local sample="$1"

    log "   -> Generating full sample list with query"
    local intermediate_list
    intermediate_list=$(mktemp)
    cat "${SAMPLE_LIST}" > "${intermediate_list}"
    echo "${sample}" >> "${intermediate_list}"

    log "   -> Generating the BED file for the first sample"
    local intermediate_bed
    intermediate_bed=$(mktemp)
    awk -F'\t' -v sample="${sample}" '$4 == sample' "${BED}" > "${intermediate_bed}"

    log "   -> Generating output VCF"
    bcftools view -Oz -o "${OUTPUT_DIR}/${sample}.vcf.gz" --samples-file "${intermediate_list}" --regions-file "${intermediate_bed}" "${VCF}"
    bcftools index "${OUTPUT_DIR}/${sample}.vcf.gz"
    log "   -> Count the number of SNPs in the VCF"
    bcftools view -H "${OUTPUT_DIR}/${sample}.vcf.gz" | wc -l > "${OUTPUT_DIR}/${sample}_SNP_count.txt"
}

log "   -> Processing individual samples"
for i in "${!SAMPLES[@]}"
do
    create_sample_vcf "${SAMPLES[$i]}"
done