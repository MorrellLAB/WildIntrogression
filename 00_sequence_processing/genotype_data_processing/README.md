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

# Updated files relative to Morex_v3


