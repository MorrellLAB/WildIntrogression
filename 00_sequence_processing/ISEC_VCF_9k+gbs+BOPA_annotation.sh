#!/bin/bash

# Input and Output Directories
INPUTDIR="/Users/jacobpacheco/Desktop/DRUM_Submission"
OUTPUTDIR="/Users/jacobpacheco/Desktop/OUTPUTTED_VCFS"

# Clean and create output directory
rm -rf "$OUTPUTDIR"
mkdir -p "$OUTPUTDIR"

# Compress and index input files if necessary
echo "Ensuring all input files are compressed and indexed..."
bgzip -c "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf" > "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz"
bcftools index -f "$INPUTDIR/wbdc_318_BOPA_morex_v3.vcf.gz"
bcftools index -f "$INPUTDIR/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz"
bcftools index -f "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz"

# Perform a three-way intersection
echo -e "\nPerforming three-way intersection..."
bcftools isec -n+3 -p "$OUTPUTDIR/three_way" \
    "$INPUTDIR/wbdc_318_BOPA_morex_v3.vcf.gz" \
    "$INPUTDIR/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz" \
    "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz"

# Compress and index the intersection file
echo -e "\nCompressing and indexing the intersection file..."
bgzip -c "$OUTPUTDIR/three_way/0000.vcf" > "$OUTPUTDIR/three_way/0000.vcf.gz"
bcftools index "$OUTPUTDIR/three_way/0000.vcf.gz"

# Annotate variants in the intersection
echo -e "\nAnnotating variants in the intersection..."
bcftools annotate \
    -a "$INPUTDIR/wbdc_318_BOPA_morex_v3.vcf.gz" \
    -c CHROM,POS,ID \
    -o "$OUTPUTDIR/temp_annotated_bopa.vcf" \
    "$OUTPUTDIR/three_way/0000.vcf.gz"

bgzip -c "$OUTPUTDIR/temp_annotated_bopa.vcf" > "$OUTPUTDIR/temp_annotated_bopa.vcf.gz"
bcftools index "$OUTPUTDIR/temp_annotated_bopa.vcf.gz"

echo -e "\nAdding 9k IDs to the intersection annotations..."
bcftools annotate \
    -a "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz" \
    -c CHROM,POS,ID \
    -o "$OUTPUTDIR/temp_annotated_final.vcf" \
    "$OUTPUTDIR/temp_annotated_bopa.vcf.gz"

bgzip -c "$OUTPUTDIR/temp_annotated_final.vcf" > "$OUTPUTDIR/temp_annotated_final.vcf.gz"
bcftools index "$OUTPUTDIR/temp_annotated_final.vcf.gz"

# Merge intersection annotations back into the original GBS file
echo -e "\nMerging annotations back into the original GBS file..."
bcftools annotate \
    -a "$OUTPUTDIR/temp_annotated_final.vcf.gz" \
    -c CHROM,POS,ID \
    -k \
    -o "$OUTPUTDIR/WBDC_GBS_final_annotated.vcf" \
    "$INPUTDIR/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz"

bgzip -c "$OUTPUTDIR/WBDC_GBS_final_annotated.vcf" > "$OUTPUTDIR/WBDC_GBS_final_annotated.vcf.gz"
bcftools index "$OUTPUTDIR/WBDC_GBS_final_annotated.vcf.gz"

# Final output message
echo -e "\nDone! Final annotated file with all variants retained is: WBDC_GBS_final_annotated.vcf.gz"


# BOPA vs GBS
echo -e "\nComparing BOPA and GBS..."
bcftools isec -n=2 -p "$OUTPUTDIR/bopa_gbs" \
    "$INPUTDIR/wbdc_318_BOPA_morex_v3.vcf.gz" \
    "$INPUTDIR/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz"
echo "Number of shared variants between BOPA and GBS: $(grep -vc '^#' "$OUTPUTDIR/bopa_gbs/0001.vcf")"

# BOPA vs 9k
echo -e "\nComparing BOPA and 9k..."
bcftools isec -n=2 -p "$OUTPUTDIR/bopa_9k" \
    "$INPUTDIR/wbdc_318_BOPA_morex_v3.vcf.gz" \
    "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz"
echo "Number of shared variants between BOPA and 9k datasets: $(grep -vc '^#' "$OUTPUTDIR/bopa_9k/0001.vcf")"

# GBS vs 9k
echo -e "\nComparing GBS and 9k..."
bcftools isec -n=2 -p "$OUTPUTDIR/gbs_9k" \
    "$INPUTDIR/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz" \
    "$INPUTDIR/9k_idt90_noRescuedSNPs.vcf.gz"
echo "Number of shared variants between GBS and 9k datasets: $(grep -vc '^#' "$OUTPUTDIR/gbs_9k/0001.vcf")"

