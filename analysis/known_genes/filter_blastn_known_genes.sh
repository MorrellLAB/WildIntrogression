#!/bin/bash

# Parse IPK Galaxy blastn results (the extended tabular format)

# User provided input arguments
FILE=~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes/blastn_domestication_genes_gene_sequence.fasta.txt_vs_morex_v3.txt

OUT_DIR=~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes

PERCENT_IDT="95"
E_VALUE="0.0001"
BIT_SCORE="60"

#-------------------
out_prefix=$(basename ${FILE} .txt)

# First pass filtering
awk -v percent_idt=${PERCENT_IDT} -v e_value=${E_VALUE} -v bit_score=${BIT_SCORE} '{ if(($3 >= percent_idt) && ($11 <= e_value) && ($12 >= bit_score)) { print } }' ${FILE} > ${OUT_DIR}/${out_prefix}.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt

# Get uniq gene ids
id_arr=($(awk '{ print $1 }' ${OUT_DIR}/${out_prefix}.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt | sort -V | uniq))
# Check count
echo ${#id_arr[@]}

if [[ -f ${OUT_DIR}/${out_prefix}.top_hit.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt ]]; then
    echo "Removing existing top hits file to start from clean slate since we are appending."
    rm ${OUT_DIR}/${out_prefix}.top_hit.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt
fi
# For each gene id, sort by highest percent identity, lowest e-value
# Then pick top "hit"
for i in "${id_arr[@]}"
do
    echo $i
    grep "$i" ${OUT_DIR}/${out_prefix}.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt | sort -k3,3nr -k11,11g | sed -n '1,1p' >> ${OUT_DIR}/${out_prefix}.top_hit.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt
done

# Convert blast tabular to BED format
# BLAST is 1-based so subtract 1 from qstart and keep qend 1-based makes it compatible with bed open interval style
# If qstart > qend (on the reverese strand), swap the two
#awk -v OFS="\t" '{ print $2,$9-1,$10,$1 }' ${OUT_DIR}/${out_prefix}.top_hit.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt | sed -e 's,lcl|,,' | tr '_' '\t' | cut -f 1-4 | sort -k1,1 -k2,2n | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' > ${OUT_DIR}/${out_prefix}.top_hits.bed
awk -v OFS="\t" '$9 > $10 { print $2,$10-1,$9,$1 } $9 < $10 { print $2,$9-1,$10,$1 }' ${OUT_DIR}/${out_prefix}.top_hit.pidt${PERCENT_IDT}_eval${E_VALUE}_bitscore${BIT_SCORE}.txt | sed -e 's,lcl|,,' | tr '_' '\t' | cut -f 1-4 | sort -k1,1 -k2,2n | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' > ${OUT_DIR}/${out_prefix}.top_hits.bed
