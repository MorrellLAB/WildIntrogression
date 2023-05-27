#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48gb
#SBATCH --tmp=22gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script converst VCF to Annovar's input file format

# Dependencies
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools

# User provided input arguments
annovar_input="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar/dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt"
# Full filepath to output directory
#   This directory should also contain the database fasta file
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar"
# This must match the *refGene.txt file
#   Example: HV_Morex_v2_HC_refGene.txt should have build version HV_Morex_v2_HC
build_version="hv_morex_v3_hc"

#-----------------
# Check if out directory exists, if not make it
mkdir -p ${out_dir}

cd ${out_dir}

annotate_variation.pl --geneanno --dbtype refGene --buildver ${build_version} ${annovar_input} ${out_dir}/
