#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=22gb
#SBATCH --tmp=12gb
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# Calculate Ts/Tv

# Dependencies
module load parallel/20210822
module load vcftools_ML/0.1.16
module load java/openjdk-17.0.2
SNPSIFT_JAR="/panfs/jay/groups/9/morrellp/shared/Software/snpEff/SnpSift.jar"

# User provided input arguments
# Prior to filtering on various annotations
#VCF1="/scratch.global/liux1299/temp_introgression/dom_and_wild_snps_polymorphic.vcf.gz"
# After filtering on various annotations
VCF2="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.vcf.gz"
VCF3="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.vcf.gz"

# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/vcf_tstv"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR}

function summarize_tstv() {
    local vcf="$1"
    local scratch_dir="$2"
    if [[ ${vcf} == *".gz"* ]]; then
        out_prefix=$(basename ${vcf} .vcf.gz)
        # Use gzip flag
        vcftools --gzvcf ${vcf} --FILTER-summary --out ${scratch_dir}/${out_prefix}
        vcftools --gzvcf ${vcf} --TsTv-summary --out ${scratch_dir}/${out_prefix}
    else
        out_prefix=$(basename ${vcf} .vcf)
        vcftools --vcf ${vcf} --FILTER-summary --out ${scratch_dir}/${out_prefix}
        vcftools --vcf ${vcf} --TsTv-summary --out ${scratch_dir}/${out_prefix}
    fi
}

export -f summarize_tstv

# Prepare output filenames
#vcf1_prefix=$(basename ${VCF1} .vcf.gz)
vcf2_prefix=$(basename ${VCF2} .vcf.gz)
vcf3_prefix=$(basename ${VCF3} .vcf.gz)

# Calculate per sample Ts/Tv
# Per sample Ts/Tv
# Redirect only summary (stdout) to output file
# VCF3 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF3} 1> ${OUT_DIR}/${vcf3_prefix}.snpsiftTsTv.txt
# VCF2 per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF2} 1> ${OUT_DIR}/${vcf2_prefix}.snpsiftTsTv.txt
# VCF1 per sample
#java -jar ${SNPSIFT_JAR} tstv ${VCF1} 1> ${OUT_DIR}/${vcf1_prefix}.snpsiftTsTv.txt

# Calculate Ts/Tv overall summary
#parallel --verbose summarize_tstv {} ${OUT_DIR} ::: ${VCF1} ${VCF2} ${VCF3}
parallel --verbose summarize_tstv {} ${OUT_DIR} ::: ${VCF2} ${VCF3}
