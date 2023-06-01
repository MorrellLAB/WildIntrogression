#!/bin/bash

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
module load bcftools/1.10.2
module load htslib/1.9
module load bedops_ML/2.4.38
module load bedtools/2.29.2

GENOMICS_GENERAL="/panfs/jay/groups/9/morrellp/liux1299/Software/genomics_general"

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz"
VCF_OUTGROUP="/panfs/jay/groups/9/morrellp/shared/Datasets/Outgroups/morex_v3_outgroups_vcf/murinum_snps_final_pseudo.vcf.gz"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/abbababa/nsgc_and_wbdc_geno"
OUT_PREFIX="dom-wild-Hmurinum_morex_v3"

WIN_SIZE="20000"
STEP_SIZE="250"

COMPLETE_SAMP_INFO="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/ABBA-BABA/Sample_Info_Complete.csv"
COMPLETE_SAMP_INFO2="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/analysis/ABBA-BABA/Sample_Info_Complete.flare_gte10percent.csv"

#-------------------
# Generate a list of regions of the VCF containing samples
# Outgroup VCF likely contains many sites that are not present in our samples
# Here we only care about sitez present in our samples
vcf_prefix=$(basename ${VCF} .vcf.gz)
zcat ${VCF} | vcf2bed | cut -f 1-4 > ${OUT_DIR}/${vcf_prefix}.bed
# Pull out sites where outgroup intersects with samples
vcf_outgroup_prefix=$(basename ${VCF_OUTGROUP} .vcf.gz)
bedtools intersect -header -wa -a ${VCF_OUTGROUP} -b ${OUT_DIR}/${vcf_prefix}.bed > ${OUT_DIR}/${vcf_outgroup_prefix}.intersect.vcf
bgzip ${OUT_DIR}/${vcf_outgroup_prefix}.intersect.vcf
tabix -p vcf --csi ${OUT_DIR}/${vcf_outgroup_prefix}.intersect.vcf.gz

# Merge vcf containing samples with vcf containing outgroup
bcftools merge -m snps -O z -o ${OUT_DIR}/${OUT_PREFIX}.vcf.gz ${OUT_DIR}/${vcf_outgroup_prefix}.intersect.vcf.gz ${VCF}

# Convert to ABBA BABA .geno.gz format
python ${GENOMICS_GENERAL}/VCF_processing/parseVCF.py -i ${OUT_DIR}/${OUT_PREFIX}.vcf.gz | bgzip > ${OUT_DIR}/${OUT_PREFIX}.geno.gz

# Follow the tutorial: https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_whole_genome/README.md
# Generate the pops.txt file with 2 columns
#    1) sample names
#    2) population name
bcftools query -l ${OUT_DIR}/${OUT_PREFIX}.vcf.gz > ${OUT_DIR}/vcf_sample_names_dom_wild_Hmurinum.txt
~/GitHub/WildIntrogression/analysis/ABBA-BABA/generate_pops_file.py ${COMPLETE_SAMP_INFO} ${OUT_DIR}/vcf_sample_names_dom_wild_Hmurinum.txt > ${OUT_DIR}/dom-wild-Hmurinum.pops.txt

# Genome-wide allele frequencies
python ${GENOMICS_GENERAL}/freq.py -g ${OUT_DIR}/${OUT_PREFIX}.geno.gz \
-p wild -p wild_introgressed -p domesticated -p outgroup \
--popsFile ${OUT_DIR}/dom-wild-Hmurinum.pops.txt --target derived \
--threads 6 \
-o ${OUT_DIR}/wild_and_dom.Hmurinum.derFreq.tsv.gz

#--------------------
# Use higher proportion introgressed SNPs threshold (10%)
# Generate the pops.txt file with 2 columns
#    1) sample names
#    2) population name
~/GitHub/WildIntrogression/analysis/ABBA-BABA/generate_pops_file.py ${COMPLETE_SAMP_INFO2} ${OUT_DIR}/vcf_sample_names_dom_wild_Hmurinum.txt > ${OUT_DIR}/dom-wild-Hmurinum.flare_gte10percent.pops.txt

# Genome-wide allele frequencies
python ${GENOMICS_GENERAL}/freq.py -g ${OUT_DIR}/${OUT_PREFIX}.geno.gz \
-p wild -p wild_introgressed -p domesticated -p outgroup \
--popsFile ${OUT_DIR}/dom-wild-Hmurinum.flare_gte10percent.pops.txt --target derived \
--threads 6 \
-o ${OUT_DIR}/wild_and_dom.Hmurinum.flare_gte10percent.derFreq.tsv.gz

#--------------------
# Sliding window ABBA-BABA
# Quantify introgression between P2 and P3
# Introgression between wild_introgressed and domesticated
python ${GENOMICS_GENERAL}/ABBABABAwindows.py \
-g ${OUT_DIR}/${OUT_PREFIX}.geno.gz -f phased \
-o ${OUT_DIR}/dom-wild-Hmurinum.ABBABABA_wild_wildintrog_dom_murinum.w${WIN_SIZE}m${STEP_SIZE}.csv.gz \
-P1 wild -P2 wild_introgressed -P3 domesticated -O outgroup \
--popsFile ${OUT_DIR}/dom-wild-Hmurinum.pops.txt -w ${WIN_SIZE} -m ${STEP_SIZE} --T 6
