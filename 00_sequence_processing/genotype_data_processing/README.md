# Genotype data processing

## Raw data

1. Wild BOPA data from Fang et al. 2014: [DRUM link](https://conservancy.umn.edu/handle/11299/181368)

```bash
# Stored locally in Dropbox
~/Dropbox/Projects/Wild_Introgression/Data/Genotyping_Data/Fang_et_al_2014_G3/WBDC_BOPA1.tsv
~/Dropbox/Projects/Wild_Introgression/Data/Genotyping_Data/Fang_et_al_2014_G3/WBDC_BOPA2.tsv
```

2. Wild barley 9K data from Nice et al. 2016

- The file "genotype.hmp.txt" is barley 9K iSelect genotyping from the Wild Barley NAM population [paper](http://www.genetics.org/content/203/3/1453).
    - The line listed as "M109" is Rasmusson. It was the recurrent parent in the experiment.

```bash
# Stored locally in Dropbox
~/Dropbox/Projects/Wild_Introgression/Data/Genotyping_Data/Nice_et_al_2016_wild_NAM_9K/genotype.hmp.txt
```

3. NSGC Core landraces and cultivated data from Poets et al. 2015 and Munoz-Amatriain et al. 2014.

NSGC Core genotyping data available in figshare: https://figshare.com/articles/dataset/Raw_Genotyping_Data_Barley_landraces_are_characterized_by_geographically_heterogeneous_genomic_origins/1468432

Download raw NSGC genotyping data.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/genotype_data_processing
sbatch download_data-NSGC_landrace_and_cultivated.sh
```

Didn't end up needing to re-download raw data. NSGC 9k genotype data was already processed and filtered. They are located here:

```bash
/panfs/jay/groups/9/morrellp/shared/Datasets/Genotyping/NSGC_Landraces
/panfs/jay/groups/9/morrellp/shared/Datasets/Genotyping/NSGC_Landraces/new_vcf_process
```

## Original files relative to Morex_v1

Li Lei generated the following files after some filtering of the BOPA genotype data.

I did convert all of the BOPA genotyping data into single vcf file. And I split the vcf file into two seperate ones: 318 WBDC accessions and 15 landraces

Please notice that:

- I only did filter the data with SNP calling probability of 0.95
- There are 74 or 33 SNPs are either missing reference or alternative in the alchemy data file

```bash
# In dir: /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/from_Li_Lei
# Combined 318 WBDC and 15 landrace vcf files
forced_flipped_forced_ref_flipped_test_BOPA.vcf

# Treat the heterozygotes as missing data
hetero_missing_318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf

# 318 WBDC vcf file
318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf

# 15 landraces vcf file
15landrace_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf
```

For Connor, `hetero_missing_318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf` is the correct file you can feed plink and set the --geno as "0.15", then plink will filter out all variants with missing call rates exceeding the provided value (0.15) to be removed.

---

## Update physical positions to Morex_v3

Previously generated VCF files (relative to an older version of the reference genome) for the genotype data are available in the directory: `/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1`. The relevant files are:

```bash
# Morex v1 files to update
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/domesticated_filtered_sorted.vcf
/panfs/jay/groups/9/morrellp/llei/Inversion/Alchemy_geno_forChaochih/hetero_missing_318WBDC__forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf
```

Updating Morex_v1 to Morex_v3 positions:

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
    /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/domesticated_filtered_sorted_clean_names.vcf \
    ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/bopa_idt95_noRescuedSNPs.vcf \
    ~/Projects/Introgressed/vcf/morex_v3 \
    domesticated_filtered_morex_v3_unsorted.vcf

# WBDC 318 samples
~/GitHub/WildIntrogression/00_sequence_processing/update_snp_pos.py \
    /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v1/hetero_missing_318WBDC__forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf \
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

## Merging VCF files

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

Merged VCF: `/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf`

## Check that strand orientation is correct relative to Morex v3

See: `tutorial_alchemy2vcf.md`

0 discordant SNPs, proceed with VCF.

## Check some summary metrics

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/genotype_data_processing
sbatch summarize_missingness.sh
sbatch check_het.sh
sbatch tstv_snps.sh
```

Check het/hom_alt ratio for any outliers (generated from bcftools stats output file with some custom calculations in script).

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/vcf_summary
# het/hom ratios
(head -n 1 merged_domesticated_and_wbdc_318_morex_v3.het_hom_ratios.txt && tail -n +2 merged_domesticated_and_wbdc_318_morex_v3.het_hom_ratios.txt | sort -k4,4nr) | head
# Observed heterozygosity
(head -n 1 merged_domesticated_and_wbdc_318_morex_v3.observed_heterozygosity.txt && tail -n +2 merged_domesticated_and_wbdc_318_morex_v3.observed_heterozygosity.txt | sort -k5,5nr) | head
```

Check missing.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/vcf_summary
sort -k5,5nr merged_domesticated_and_wbdc_318_morex_v3.vcf_missingness.imiss | head
sort -k6,6nr merged_domesticated_and_wbdc_318_morex_v3.vcf_missingness.lmiss | head
```

## Filter VCFs

Filter WBDC BOPA and NSGC 9K VCFs before merging (after merging causes some overly stringent filtering for SNPs present in one dataset but not in the other due to different genotyping sets).

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/genotype_data_processing
./Filter_VCF_wbdc_bopa.sh
./Filter_VCF_nsgc.sh
```

Check summary metrics after filtering.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/genotype_data_processing
./check_het-wbdc.sh
./check_het-nsgc.sh

./summarize_missingness-wbdc.sh
./summarize_missingness-nsgc.sh

./tstv_snps-wbdc.sh
./tstv_snps-nsgc.sh
```

Check het/hom_alt ratio for any outliers (generated from bcftools stats output file with some custom calculations in script).

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/vcf_summary
# het/hom ratios
(head -n 1 domesticated_snps.polymorphic.filt_miss_het.het_hom_ratios.txt && tail -n +2 domesticated_snps.polymorphic.filt_miss_het.het_hom_ratios.txt | sort -k4,4nr) | head

(head -n 1 wbdc_bopa_snps.polymorphic.filt_miss_het.het_hom_ratios.txt && tail -n +2 wbdc_bopa_snps.polymorphic.filt_miss_het.het_hom_ratios.txt | sort -k4,4nr) | head

# Observed heterozygosity
(head -n 1 domesticated_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt && tail -n +2 domesticated_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt | sort -k5,5nr) | head

(head -n 1 wbdc_bopa_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt && tail -n +2 wbdc_bopa_snps.polymorphic.filt_miss_het.observed_heterozygosity.txt | sort -k5,5nr) | head
```

Check missing.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/vcf_summary
sort -k5,5nr domesticated_snps.polymorphic.filt_miss_het_missingness.imiss | head

sort -k6,6nr domesticated_snps.polymorphic.filt_miss_het_missingness.lmiss | head

sort -k5,5nr wbdc_bopa_snps.polymorphic.filt_miss_het_missingness.imiss | head

sort -k6,6nr wbdc_bopa_snps.polymorphic.filt_miss_het_missingness.lmiss | head
```

## Merge filtered VCF files

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
# Merge VCF files
bcftools merge -m id -O z -o nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz domesticated_snps.polymorphic.filt_miss_het.vcf.gz wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz
# Index
tabix -p vcf --csi nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz
```

Filtered VCF filepaths:

```bash
# domesticated accessions (landraces and cultivars) only
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/domesticated_snps.polymorphic.filt_miss_het.vcf.gz
# wild accessions only
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz

# domesticated and wild accessions merged into one file
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz
```
