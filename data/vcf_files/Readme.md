I did convert all of the BOPA genotyping data into single vcf file. And I split the vcf file into two seperate ones: 318 WBDC accessions and 15 landraces

Please notice that:

- I only did filter the data with SNP calling probability of 0.95

- There are 74 or 33 SNPs are either missing reference or alternative in the alchemy data file


`forced_flipped_forced_ref_flipped_test_BOPA.vcf`: combined 318 WBDC and 15 landraces vcf file

`hetero_missing_318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf`: treat the heterozygotes as missing data

`318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf`: 318 WBDC vcf file

`15landrace_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf`: 15 landraces vcf file


For Connor, `hetero_missing_318WBDC_forced_flipped_forced_ref_flipped_test_BOPA.recode.vcf` is the correct file you can feed plink and set the --geno as "0.15", then plink will filter out all variants with missing call rates exceeding the provided value (0.15) to be removed.
