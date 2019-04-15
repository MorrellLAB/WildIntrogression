Connor Depies, 29/05/2018, St. Paul Minnesota
ibs, sequence_handling are not mine. ibs has a readme inside it.

hmp_to_PLINK.py File Tom Kono -Dependencies sys, pprint
"""Convert from the HMP format from T3 into PLINK PED files. Will order markers
according to a supplied PLINK Map file, and check alleles for RC errors with
a VCF. Takes three agruments:
    1) T3 hmp file
    2) PLINK MAP
    3) BOPA VCF"""
merge_vcf.pl

##by Li Lei, 20170809, in St. Paul;
#this is to merge vcf files for Connor. Since vcf-merge can not handle the the position out of 500 MB. 
# The prerequest is to make sure all of the SNPs in the 1-9 columns are identical. 
You can run this commandline to make sure if they are the same:  
diff -y <(grep -v "#" sorted_filtered_round2.vcf|cut -f 1,2,3,4,5,9) <(grep -v "#" sorted_filtered_NAM_9k.vcf|cut -f 1,2,3,4,5,9)|less -S
#usage: perl /panfs/roc/groups/9/morrellp/llei/Introgressed_line/script/merge_vcf.pl sorted_filtered_round2.vcf sorted_filtered_NAM_9k.vcf >merged_NAM_9k_round2.vcf

Matrix_Manipulation Directory -Chaochih, Connor: R-codes for working with matrices, used for neighbor joining trees
split_by_snips.sh

#Connor Depies
#10/30/2017
#Program to split a vcf file into 100 snp intervals sliding by 25 snps and analyze their IBD with plink
#1) sorted VCF file
Dependencies plink/1.90b

The following 2 files are more notes than scripts. You can run them, but they only really work on the files which I used them on.

mergehmpvcf.sh
#Connor Depies Aug 17, 2017
#This program is very specific to the data I have been working with, 
so I do not suggest anyone use it for anything other than the following files in order: 
iSelect_9k.map, sorted_all_9k_masked_90idt.vcf, genome.hmp.txt and NAM_9k_Final.vcf
#This program takes 4 inputs:
#1 ) A map file
#2 ) A partial vcf file with the same SNPs as the 4
#3 ) A genotype file
#These three files are filtered of unshared snps, then converted to a full vcf file. 
Then the program fixes flip and force errors (see Li's tutorial on the matter), 
removes snps which don't have both alleles.
#4 ) A full vcf file
# The paths to two programs must also be inputted: one
 of them written by Tom to convert hmp to plink, the other by Li to merge vcf files:
#We then merge the vcf files
Dependencies plink/1.90b python3 vcftools_ML/0.1.14 samtools_ML/1.3.1 htslib_ML/1.4.0 perl/5.14.2 
Script dependencies hmp_to_PLINK.py merge_vcf.pl

flips_and_forces.sh
Basically takes the latter part of mergehmpvcf.sh 
#Connor Depies
#08/03/2018
# Fixes flip and force errors so that two vcf files which have already been sorted and filtered so that they have the same SNPs 
# will also have the same major and minor allele assignments and no missing values.
# Takes two sorted vcf files.
# Takes two strings for names of of output files
# Takes one more string for directory in which to place intermediate files
# Takes the number of lines of each of the two headers
# Adapted from Li Lei's Tutorial
Dependencies plink/1.90b
vcftools_ML/0.1.14

