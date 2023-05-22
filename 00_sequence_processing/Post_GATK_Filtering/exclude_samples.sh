#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --tmp=10gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2

VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.vcf.gz"
SAMPLES_TO_EXCLUDE="/panfs/jay/groups/9/morrellp/liux1299/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering/excluded_samples.txt"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered"

out_prefix=$(basename ${VCF} .vcf.gz)

# In dir: ~/Alignments/introgression_project/all_dom_and_wild/Filtered
bcftools view --samples-file ^${SAMPLES_TO_EXCLUDE} ${VCF} -O z -o ${OUT_DIR}/${out_prefix}.final.vcf.gz
