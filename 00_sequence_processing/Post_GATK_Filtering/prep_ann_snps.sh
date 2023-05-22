#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --tmp=8gb
#SBATCH -t 80:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script prepares input files for scikit_allel_plot_prep.sh and variant_analysis
# Adds known vs novel annotations to VCF
# Adds filtered vs retained_filt1 vs retained_filt2 annotations to VCF

# Dependencies
module load bedops_ML/2.4.38
module load bedtools/2.29.2
module load bcftools/1.10.2
module load htslib/1.9
module load python3/3.8.3_anaconda2020.07_mamba
export PATH=${PATH}:"/panfs/jay/groups/9/morrellp/liux1299/GitHub/Barley_Mutated/01_snp_filtering/Post_GATK_Filtering"

# User provided input arguments
VCF_UNFILT="/scratch.global/liux1299/temp_introgression/dom_and_wild_snps_polymorphic.vcf.gz"
FILT1_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.vcf.gz"
FILT2_VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.vcf.gz"
# VCFs of known variants
VCF_SANGER="/panfs/jay/groups/9/morrellp/pmorrell/Workshop/Selective_Sweeps/Sanger/Morex_v3_processed/17_barley_sanger_loci_Morex_v3_parts.vcf"
VCF_BOPA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/bopa_idt95_noRescuedSNPs_partsRef.vcf"
VCF_9K="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/9k_idt95_noRescuedSNPs_partsRef.vcf"
VCF_50K="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/bopa_9k_50k/50k_idt95_noRescuedSNPs_partsRef.vcf"
# Output file prefix
KNOWN_OUT_PREFIX="sanger_bopa_9k_50k_idt95_noRescuedSNPs_partsRef"
VCF_OUT_PREFIX="dom_and_wild_snps"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/scikit_allel_files"
SCRATCH_DIR="/scratch.global/liux1299/temp_introgression/scikit_allel_files"
# Select chr1H only
CHROM_PREFIX="chr1H"

#-----------------------
mkdir -p ${OUT_DIR} ${SCRATCH_DIR}

function vcf_to_bed() {
    local vcf="$1"
    local out_dir="$2"
    if [[ "${vcf}" == *"gz"* ]]; then
        # We are working with gzipped VCF file
        out_prefix=$(basename ${vcf} .vcf.gz)
    else
        # We are working with uncompressed vcf
        out_prefix=$(basename ${vcf} .vcf)
    fi
    # Convert VCF to bed format
    vcf2bed < ${vcf} > ${out_dir}/${out_prefix}.bed
}

export -f vcf_to_bed

# Get the basename of files
bn_sanger=$(basename ${VCF_SANGER} .vcf)
bn_bopa=$(basename ${VCF_BOPA} .vcf)
bn_9k=$(basename ${VCF_9K} .vcf)
bn_50k=$(basename ${VCF_50K} .vcf)

# Convert VCF to bed format for all known variants
vcf_to_bed ${VCF_BOPA} ${OUT_DIR}
vcf_to_bed ${VCF_9K} ${OUT_DIR}
vcf_to_bed ${VCF_50K} ${OUT_DIR}
# For Sanger, we only want SNPs so we'll process this one differently from rest
vcf2bed < ${VCF_SANGER} | grep "SNP" > ${OUT_DIR}/${bn_sanger}-snps.bed

# Concatenate BED of known variants and get unique, sorted positions
cat ${OUT_DIR}/${bn_sanger}-snps.bed ${OUT_DIR}/${bn_bopa}.bed ${OUT_DIR}/${bn_9k}.bed ${OUT_DIR}/${bn_50k}.bed | cut -f 1-3 | sort -k1,1 -k2,2n | uniq > ${OUT_DIR}/${KNOWN_OUT_PREFIX}.bed

bedtools intersect -header -wa -a ${VCF_UNFILT} -b ${OUT_DIR}/${KNOWN_OUT_PREFIX}.bed > ${OUT_DIR}/known_x_${VCF_OUT_PREFIX}.vcf

# Create new tab delimited file with new annotation added
# Add known vs novel variant annotations
vcf_to_bcftools_ann_tsv-known.py ${VCF_UNFILT} ${OUT_DIR}/known_x_${VCF_OUT_PREFIX}.vcf > ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.txt

# Next, add filtered vs retained_filt1 vs retained_filt2 annotations
vcf_to_bcftools_ann_tsv-filtered.py ${VCF_UNFILT} ${FILT1_VCF} ${FILT2_VCF} > ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.txt

# bgzip and tabix index before running bcftools annotate
bgzip --force ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.txt
tabix -p vcf ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.txt.gz
bgzip --force ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.txt
tabix -p vcf ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.txt.gz

# Prepare .hdr files containing VCF header lines to be added to output vcf
echo '##INFO=<ID=VAR_KNOWN,Number=1,Type=String,Description="Custom variant category, novel vs known">' > ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.hdr
echo '##INFO=<ID=VAR_FILT,Number=1,Type=String,Description="Custom variant category corresponding to each level of filtering, filtered vs retained_filt1 vs retained_filt2">' > ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.hdr

# Add new TAG for known vs novel variant to INFO field
bcftools annotate -a ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.txt.gz -h ${OUT_DIR}/ann_known-${VCF_OUT_PREFIX}.hdr -c CHROM,POS,REF,ALT,-,VAR_KNOWN ${VCF_UNFILT} -O z -o ${SCRATCH_DIR}/ann_known-${VCF_OUT_PREFIX}.vcf.gz

# Add filtered TAG for filtered vs retained_filt1 vs retained_filt2 to INFO field
bcftools annotate -a ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.txt.gz -h ${OUT_DIR}/ann_filt-${VCF_OUT_PREFIX}.hdr -c CHROM,POS,REF,ALT,-,VAR_FILT ${SCRATCH_DIR}/ann_known-${VCF_OUT_PREFIX}.vcf.gz -O z -o ${OUT_DIR}/ann_known_and_filt-${VCF_OUT_PREFIX}.vcf.gz

# Select chr1H for visual exploration
bcftools view -t "chr1H_part1,chr1H_part2" ${OUT_DIR}/ann_known_and_filt-${VCF_OUT_PREFIX}.vcf.gz -O z -o ${OUT_DIR}/ann_known_and_filt-${CHROM_PREFIX}_${VCF_OUT_PREFIX}.vcf.gz

# Pull out known variants that got filtered out
bcftools filter -i 'VAR_KNOWN="known" & VAR_FILT!="retained_filt2"' ${OUT_DIR}/ann_known_and_filt-${VCF_OUT_PREFIX}.vcf.gz -O z -o ${OUT_DIR}/only_known-filtered-retained_filt1-${VCF_OUT_PREFIX}.vcf.gz

bcftools filter -i 'VAR_KNOWN="known"' ${OUT_DIR}/ann_known_and_filt-${VCF_OUT_PREFIX}.vcf.gz -O z -o ${OUT_DIR}/only_known-${VCF_OUT_PREFIX}.vcf.gz
