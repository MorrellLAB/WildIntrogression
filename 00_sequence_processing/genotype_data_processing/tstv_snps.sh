#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=8gb
#SBATCH -t 00:10:00
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
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf"
# Output directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/vcf_tstv"

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
vcf_prefix=$(basename ${VCF} .vcf.gz)

# Calculate per sample Ts/Tv
# Per sample Ts/Tv
# Redirect only summary (stdout) to output file
# VCF per sample
java -jar ${SNPSIFT_JAR} tstv ${VCF} 1> ${OUT_DIR}/${vcf_prefix}.snpsiftTsTv.txt

# Calculate Ts/Tv overall summary
summarize_tstv ${VCF} ${OUT_DIR}
