#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12gb
#SBATCH --tmp=6gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Converts Annovar output file to BAD_Mutations subs files

# Dependencies
module load python3/3.8.3_anaconda2020.07_mamba
# Directory containing ANNOVAR_to_subs.py script
export PATH=${PATH}:~/GitHub/WildIntrogression/analysis/Annovar

# User provided input arguments
ANNOVAR_EXONIC_VAR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar/dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt.exonic_variant_function"
ANNOVAR_ALL_EFFECTS="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar/dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt.variant_function"
SUBDIR_BN="dom_and_wild_snps"
OUT_NAME="dom_and_wild_snps"
# Full path to out directory
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations/annovar_to_subs"

SCRIPT_ANNOVAR_TO_EFFECTS=~/GitHub/Barley_Mutated/02_analysis/Annovar/ANNOVAR_To_Effects.py

#-----------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/${SUBDIR_BN}

# Get directory of annovar file
annovar_outputs_dir=$(dirname ${ANNOVAR_EXONIC_VAR})
out_prefix=$(basename ${ANNOVAR_EXONIC_VAR})

# Pull nonsynonymous variants as defined by Sequence Ontology (http://www.sequenceontology.org/)
# nonsynonymous SNV in Annovar's report is a missense variant
# See table in: https://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-2-refseq-gene-annotation
grep 'nonsynonymous SNV\|stopgain\|stoploss' ${ANNOVAR_EXONIC_VAR} | uniq > ${annovar_outputs_dir}/${out_prefix}.nonsyn

# Annovar .exonic_variant_function file format
# Column 1: line # in original input file
# Column 2: functional consequences
# Column 3: gene name : transcript identifier : sequence change in corresponding transcript
#   For nomenclature, see: https://varnomen.hgvs.org/bg-material/simple/
ANNOVAR_to_subs.py ${annovar_outputs_dir}/${out_prefix}.nonsyn "${OUT_DIR}/${SUBDIR_BN}/${OUT_NAME}.long_subs.txt" "${OUT_DIR}/${SUBDIR_BN}"
