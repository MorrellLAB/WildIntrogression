I did convert all of the BOPA genotyping data into single vcf file. And I split the vcf file into two seperate ones: 318 WBDC accessions and 15 landraces

Please notice that:

- I only did filer the data with SNP calling probability of 0.95

- There are 74 or 33 SNPs are either missing reference or alternative in the alchemy data file

- You may need to treat the heterozygotes as missing and filter out the SNPs with certain missing individuals.
