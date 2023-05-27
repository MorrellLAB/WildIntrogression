#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=6gb
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script convert VCF to Annovar's input file format

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
vcf_file="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.final.vcf.gz"
# Full filepath to output directory
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar"

#-----------------
# Generate output file pfrefix
out_prefix=$(basename ${vcf_file} .vcf.gz)

convert2annovar.pl --format vcf4 --allsample --withfreq --includeinfo --outfile ${out_dir}/${out_prefix}_annovar_input.txt ${vcf_file}
