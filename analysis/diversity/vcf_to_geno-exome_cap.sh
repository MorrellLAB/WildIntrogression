#!/bin/bash

module load python3/3.8.3_anaconda2020.07_mamba
module load htslib/1.9
module load bcftools/1.10.2
GGDIR="/panfs/jay/groups/9/morrellp/liux1299/Software/genomics_general"

VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.final.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/diversity"

#------------
mkdir -p ${OUT_DIR}

out_prefix=$(basename $VCF .vcf.gz)

# Convert VCF to .geno.gz format
python ${GGDIR}/VCF_processing/parseVCF.py -i ${VCF} | bgzip > ${OUT_DIR}/${out_prefix}.geno.gz
