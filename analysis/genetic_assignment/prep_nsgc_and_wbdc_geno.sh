#!/bin/bash

# Dependencies
module load bcftools/1.10.2
module load plink/1.90b6.10
module load python3/3.8.3_anaconda2020.07_mamba

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data_morex_v3"
OUT_PREFIX="nsgc_and_wbdc_geno_snps_morex_v3"

POP_INFO="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/genetic_assignment/pop_info_nsgc_and_wild_geno.txt"

GENETIC_MAP_BOPA="/panfs/jay/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt"
cM_COL_NUM_BOPA="3"
ID_COL_NUM_BOPA="2"
GENETIC_MAP_9k="/panfs/jay/groups/9/morrellp/liux1299/GitHub/morex_reference/genetic_maps/9k_iSelect_Genetic_Map.txt"
cM_COL_NUM_9k="4"
ID_COL_NUM_9k="1"

# Parameters for Plink LD pruning
WIN_SIZE="50"
STEP_SIZE="10"
R2_THRESHOLD="0.1"

#----------
# Rename samples in VCF to include population info
bcftools reheader -s ${POP_INFO} -o ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.vcf.gz ${VCF}

# Convert VCF to Plink 1.9 MAP/PED/BED/BIM files
# Add cM positions to Plink .map file
plink --vcf ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.vcf.gz --keep-allele-order --allow-extra-chr --id-delim ":" --make-bed --update-cm ${GENETIC_MAP_BOPA} ${cM_COL_NUM_BOPA} ${ID_COL_NUM_BOPA} --recode --out ${OUT_DIR}/${OUT_PREFIX}_wPopInfo
# Second round of updating cM positions in Plink .map file with 9k SNPs
plink --bfile ${OUT_DIR}/${OUT_PREFIX}_wPopInfo -keep-allele-order --allow-extra-chr --update-cm ${GENETIC_MAP_9k} ${cM_COL_NUM_9k} ${ID_COL_NUM_9k} --recode --out ${OUT_DIR}/${OUT_PREFIX}_wPopInfo

# In the MAP file, make sure chromosomes are designated with their chromosome number only and don't have additional characters. Example, `chr1H` gets changed to `1`.
# We'll use sed to do an in-place edit (so we don't have to keep renaming the output files to match plink's file set naming system)
for i in $(seq 1 7)
do
    # Replace for MAP and BIM files
    sed -i "s/chr${i}H/${i}/g" ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.map
    sed -i "s/chr${i}H/${i}/g" ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.bim
done

# LD pruning with Plink
plink --bfile ${OUT_DIR}/${OUT_PREFIX}_wPopInfo --indep-pairwise ${WIN_SIZE} ${STEP_SIZE} ${R2_THRESHOLD} --out ${OUT_DIR}/${OUT_PREFIX}_wPopInfo
# 977 of 2619 variants removed, 1,642 variants remaining for analysis

# Generate BED/BIM files after LD pruning
plink --bfile ${OUT_DIR}/${OUT_PREFIX}_wPopInfo --extract ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.prune.in --make-bed --out ${OUT_DIR}/${OUT_PREFIX}_wPopInfo.pruned
