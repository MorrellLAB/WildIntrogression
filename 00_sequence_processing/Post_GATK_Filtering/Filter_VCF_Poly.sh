#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=9gb
#SBATCH --tmp=8gb
#SBATCH -t 60:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Dependencies
module load bcftools/1.10.2
module load gatk/4.1.2
module load htslib/1.9

# User provided input arguments
VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Variant_Recalibrator/dom_and_wild_snps.recalibrated.pass_sites.vcf.gz"
OUT_PREFIX="dom_and_wild_snps"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered"
# Scratch directory used to store intermediate files that don't need to be kept long term
SCRATCH_DIR="/scratch.global/liux1299/temp_introgression"

#-----------------
# Check the out dir and scratch dir exist, if not make them
mkdir -p ${OUT_DIR} ${SCRATCH_DIR}

# Check the number of sites in starting VCF
num_sites_vcf=$(zgrep -v "#" ${VCF} | wc -l)
# Append the number of sites remaining to file
printf "${VCF}\t${num_sites_vcf}\n" >> ${OUT_DIR}/${OUT_PREFIX}_num_sites.log

# Remove sites that aren't polymorphic (minor allele count of 0) and unused alternate alleles.
echo "Removing sites that aren't polymorphic and unused alternate alleles..."
gatk SelectVariants \
    -V "${VCF}" \
    --exclude-non-variants true \
    --remove-unused-alternates true \
    -O "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz" \
    --tmp-dir ${SCRATCH_DIR}
echo "Done removing sites that aren't polymorphic and unused alternate alleles."
# Get the number of sites left after filtering
num_sites_poly=$(zgrep -v "#" ${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz | wc -l)
# Append the number of sites remaining to file
printf "${SCRATCH_DIR}/${OUT_PREFIX}_polymorphic.vcf.gz\t${num_sites_poly}\n" >> ${OUT_DIR}/${OUT_PREFIX}_num_sites.log
