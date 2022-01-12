# Sequence Processing

Sequence alignment and SNP calling of the WBDC samples were performed as part of another larger project. The configs are available in the GitHub repository: https://github.com/MorrellLAB/Barley_Inversions/tree/master/00_sequence_processing/configs_morex_v3.

This directory also includes other scripts related to data processing steps (e.g., updating physical positions of VCF files, etc.).

---

### Update Physcial Positions

Previously generated VCF files (relative to an older version of the reference genome) for the genotype data are available in the directory: `/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1`

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v1
# First, cleanup sample names (so they are not "doubled" names)
# Make sure to remove trailing tab
grep "#CHROM" domesticated_filtered_sorted.vcf | tr '\t' '\n' | cut -d'_' -f 1 | tr '\n' '\t' | sed -e '$a\' | sed 's/[\t]*$//' > temp_domesticated_fixed_sample_names.txt
# Confirm samples are in the same order
grep "#CHROM" domesticated_filtered_sorted.vcf | tr '\t' '\n' | cut -d'_' -f 1 > temp_new_dom_samp_names.txt
grep "#CHROM" domesticated_filtered_sorted.vcf | tr '\t' '\n' > temp_original_dom_samp_names.txt
diff -y temp_original_dom_samp_names.txt temp_new_dom_samp_names.txt
# Generate new file with the fixed sample names
grep "##" domesticated_filtered_sorted.vcf > domesticated_filtered_sorted_clean_names.vcf
cat temp_domesticated_fixed_sample_names.txt >> domesticated_filtered_sorted_clean_names.vcf
grep -v "#" domesticated_filtered_sorted.vcf >> domesticated_filtered_sorted_clean_names.vcf
# Check column counts one more time
grep "#CHROM" domesticated_filtered_sorted_clean_names.vcf | tr '\t' '\n' | wc -l
grep -v "#" domesticated_filtered_sorted_clean_names.vcf | head -n 1 | tr '\t' '\n' | wc -l

# Update positions
module load python3/3.8.3_anaconda2020.07_mamba
# Domesticated samples
~/GitHub/WildIntrogression/00_sequence_processing/update_snp_pos.py \
    /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/domesticated_filtered_sorted_clean_names.vcf \
    ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/bopa_idt95_noRescuedSNPs.vcf \
    ~/Projects/Introgressed/vcf/morex_v3 \
    domesticated_filtered_morex_v3_unsorted.vcf

# WBDC 318 samples
~/GitHub/WildIntrogression/00_sequence_processing/update_snp_pos.py \
    /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/hetero_missing_318WBDC__forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf \
    ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/bopa_idt95_noRescuedSNPs.vcf \
    ~/Projects/Introgressed/vcf/morex_v3 \
    wbdc_318_BOPA_morex_v3_unsorted.vcf
```

Sort the VCF files.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
module load bcftools/1.10.2
# Domesticated samples
bcftools sort domesticated_filtered_morex_v3_unsorted.vcf -O v -o domesticated_filtered_morex_v3.vcf

# WBDC 318 samples
bcftools sort wbdc_318_BOPA_morex_v3_unsorted.vcf -O v -o wbdc_318_BOPA_morex_v3.vcf

# Clean up intermediate files
rm domesticated_filtered_morex_v3_unsorted.vcf
rm wbdc_318_BOPA_morex_v3_unsorted.vcf

# Check the number of SNPs where we don't have the Morex v3 positions
# due to difficulty with the BOPA SNP alignment using SNP_Utils
grep -v "#" missing_snps_list_domesticated_filtered_morex_v3_unsorted.vcf | wc -l
    43
grep -v "#" missing_snps_list_wbdc_318_BOPA_morex_v3_unsorted.vcf | wc -l
    75
```

### Merging VCF files

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
# Dependencies
module load htslib/1.9
module load bcftools/1.10.2
# bcftools requires bgzipped vcf files for merge
bgzip -c domesticated_filtered_morex_v3.vcf > domesticated_filtered_morex_v3.vcf.gz
tabix -p vcf --csi domesticated_filtered_morex_v3.vcf.gz
bgzip -c wbdc_318_BOPA_morex_v3.vcf > wbdc_318_BOPA_morex_v3.vcf.gz
tabix -p vcf --csi wbdc_318_BOPA_morex_v3.vcf.gz

# Merge VCF files
bcftools merge -m id -O v -o merged_domesticated_and_wbdc_318_morex_v3.vcf domesticated_filtered_morex_v3.vcf.gz wbdc_318_BOPA_morex_v3.vcf.gz
# Bgzip and index
bgzip -c merged_domesticated_and_wbdc_318_morex_v3.vcf > merged_domesticated_and_wbdc_318_morex_v3.vcf.gz
tabix -p vcf --csi merged_domesticated_and_wbdc_318_morex_v3.vcf.gz
```

Merged VCF: `/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf`
