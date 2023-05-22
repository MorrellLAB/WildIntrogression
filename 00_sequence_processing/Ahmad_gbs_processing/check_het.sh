#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8gb
#SBATCH --tmp=6gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

# Dependencies
module load bcftools/1.10.2
module load texlive/20131202

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3/final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3/vcf_summary"
# Arguments for bcftools stats
REF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"

SCRIPT_BCFTOOLS_STATS=~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering/run_bcftools_stats.sh

#---------------
source ${SCRIPT_BCFTOOLS_STATS}

# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR}

# Generate bcftools stats
generate_stats ${VCF} ${REF} ${OUT_DIR}
