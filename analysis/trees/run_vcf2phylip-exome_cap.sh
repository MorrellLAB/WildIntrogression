#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --tmp=2gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load python3/3.8.3_anaconda2020.07_mamba
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/liux1299/Software/vcf2phylip

VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.final.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees"

#------------------
mkdir -p ${OUT_DIR}

vcf2phylip.py --input ${VCF} --output-folder ${OUT_DIR}
