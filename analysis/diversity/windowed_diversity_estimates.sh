#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --tmp=8gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load python3/3.8.3_anaconda2020.07_mamba
GGDIR="/panfs/jay/groups/9/morrellp/liux1299/Software/genomics_general"

VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/diversity"
POPS_FILE=~/GitHub/WildIntrogression/analysis/ABBA-BABA/dom-wild.pops.txt

WIN_SIZE="10000"
STEP_SIZE="250"
MIN_GS="50"

out_prefix=$(basename $VCF .vcf.gz)

THREADS=8

# Calculate sliding window diversity, we're interested in fst
python ${GGDIR}/popgenWindows.py -w ${WIN_SIZE} -s ${STEP_SIZE} -m ${MIN_GS} -g ${OUT_DIR}/${out_prefix}.geno.gz -o ${OUT_DIR}/${out_prefix}.w${WIN_SIZE}.m${STEP_SIZE}.csv.gz -f phased --addWindowID --threads ${THREADS} -p "domesticated" -p "wild" -p "wild_introgressed" --popsFile ${POPS_FILE}
