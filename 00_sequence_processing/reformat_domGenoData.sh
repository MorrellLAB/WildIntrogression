#!/bin/bash

set -e
set -o pipefail

# Log of commands run to get domesticated Genotype data VCF sample names to match sample names
#   in adapter trimmed fastq files for exome data

# We have already done some reformatting, but need to do a little more.
#   Example:
#       dom geno data VCF: CIho_497 (original: CIho497)
#       exome capture trimmed fastq sample name: CIho_00497
DOM_GENO_SAMPLE_NAMES="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/temp_reformatted_accession_names_dom_genoData.txt"
OUT_FILE="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/temp_reformatted2_accession_names_dom_genoData.txt"

#-------------------
for i in $(cat ${DOM_GENO_SAMPLE_NAMES})
do
    # Deal with sample names that begin with CIho
    if [[ ${i} == "CIho"* ]] && [[ $(echo ${i} | wc -c) == "9" ]]; then
        # add two zeros to number
        # Example: CIho_497 should changed to CIho_00497
        echo ${i} | sed 's,_,_00,' >> ${OUT_FILE}
    elif [[ ${i} == "CIho"* ]] && [[ $(echo ${i} | wc -c) == "10" ]]; then
        # add one zero to number
        # Example: CIho_1551 should be changed to CIho_01551
        echo ${i} | sed 's,_,_0,' >> ${OUT_FILE}
    # Next, deal with sample names that begin with PI
    elif [[ ${i} == "PI"* ]] && [[ $(echo ${i} | wc -c) == "8" ]]; then
        echo ${i} | sed 's,_,_00,' >> ${OUT_FILE}
    elif [[ ${i} == "PI"* ]] && [[ $(echo ${i} | wc -c) == "9" ]]; then
        echo ${i} | sed 's,_,_0,' >> ${OUT_FILE}
    else
        echo ${i} >> ${OUT_FILE}
    fi
done
