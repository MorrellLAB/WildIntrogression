#!/bin/bash

#   Chaochih Liu - May 11, 2018

#   Dependencies
#   GNU datamash

PHYS_POS_SNP_DIR=/Users/chaochih/Dropbox/Projects/Wild_Introgression/Analyses/IBS/Results/snp_windows_filtered_physPos
FILT_PI_HAT_LIST=/Users/chaochih/Dropbox/Projects/Wild_Introgression/Analyses/IBS/Results/reformatted_wild_vs_cult_introg/filtered_pi_hat_1.0/piHat_1.0_list.txt
OUT_DIR=/Users/chaochih/Dropbox/Projects/Wild_Introgression/Analyses/IBS/Results/identify_runs_of_ibs

#   Check if output directories exist, if not make them
mkdir -p ${OUT_DIR}
mkdir -p ${OUT_DIR}/temp

#   Store file list in an array
PI_HAT_ARRAY=($(cat ${FILT_PI_HAT_LIST}))

#   Create list of patterns to search
for i in ${PI_HAT_ARRAY[@]}
do
    out_prefix=$(basename ${i} .txt | cut -d'_' -f 1,2,3)
    #   First generate file of only pairs of individuals (Ind1 and Ind2 columns)
    #   Skip header line to make our lives easier
    awk '{ print $1"\t"$2 }' ${i} | tail -n +2 > ${OUT_DIR}/temp/${out_prefix}_indPairs.txt
    repeat=$(cat ${OUT_DIR}/temp/${out_prefix}_indPairs.txt | wc -l) # number of rows (pairs) we have
    #   Extract chrom and snp interval start-end info from filename
    #   (Ex filename: chr1_snpsList_0-100_reformat_genome_z2_0.90_wild_vs_cult_piHat_1.0.txt)
    chrom=$(basename ${i} .txt | cut -d'_' -f 1)
    snp_int=$(basename ${i} .txt | cut -d'_' -f 3)
    #   Print each row of chromosome and snp interval start-end n times
    for r in $(seq 1 ${repeat})
    do
        printf "${chrom}\t${snp_int}\n%.0s" >> ${OUT_DIR}/temp/tmp_${out_prefix}_chr_snpInt_metadata.txt
    done

    #   Build a new file by combining individual pairs file and metadata file columns
    paste -d'\t' ${OUT_DIR}/temp/${out_prefix}_indPairs.txt ${OUT_DIR}/temp/tmp_${out_prefix}_chr_snpInt_metadata.txt > ${OUT_DIR}/temp/${out_prefix}_chr_snpInt.txt
done



#   Go into temp directory and create file list
cd ${OUT_DIR}/temp
find $(pwd) -name "*_chr_snpInt.txt" | sort -V > tmp_chr_snpInt_list.txt
FILE_LIST_PATH=$(find $(pwd) -name tmp_chr_snpInt_list.txt)
CHR_INT_ARRAY=($(cat ${FILE_LIST_PATH}))

for i in ${CHR_INT_ARRAY[@]}
do
    shared_prefix=$(basename ${i} .txt | cut -d'_' -f 1,2,3)
    minimum=$(tail -n +2 ${PHYS_POS_SNP_DIR}/${shared_prefix}_sorted_physPos_noNA.txt | datamash min 3)
    maximum=$(tail -n +2 ${PHYS_POS_SNP_DIR}/${shared_prefix}_sorted_physPos_noNA.txt | datamash max 3)
    phys_size=$(expr ${maximum} - ${minimum})
    repeat=$(cat ${OUT_DIR}/temp/${shared_prefix}_chr_snpInt.txt | wc -l) # number of rows we have
    for r in $(seq 1 ${repeat})
    do
        printf "${minimum}\t${maximum}\t${phys_size}\n%.0s" >> ${OUT_DIR}/temp/tmp_${shared_prefix}_phys_size.txt
    done

    #   Build a new file by combining metadata columns and physical size columns
    #   Create header first
    printf "Ind1\tInd2\tChr\tSNP_Win\tPhysPos_Start\tPhysPos_End\tInt_Phys_Size\n" > ${OUT_DIR}/${shared_prefix}_win_physSize.txt
    paste -d'\t' ${OUT_DIR}/temp/${shared_prefix}_chr_snpInt.txt ${OUT_DIR}/temp/tmp_${shared_prefix}_phys_size.txt > ${OUT_DIR}/${shared_prefix}_win_physSize.txt
done

#   Concatenate all chromosomes
cd ${OUT_DIR}
printf "Ind1\tInd2\tChr\tSNP_Win\tPhysPos_Start\tPhysPos_End\tInt_Phys_Size\n" > ${OUT_DIR}/all_chr_combined_metadata.txt
cat *_win_physSize.txt | sort -V -k5,5 | sort -k2,2 | sort -k3,3 >> ${OUT_DIR}/all_chr_combined_metadata.txt
