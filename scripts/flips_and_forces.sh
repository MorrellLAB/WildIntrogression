#!/bin/env bash
#Connor Depies
#08/03/2018
# Fixes flip and force errors so that two vcf files which have already been sorted and filtered so that they have the same SNPs 
# will also have the same major and minor allele assignments and no missing values.
# Takes two vcf files.
# Takes two strings for names of of output files
# Takes one more string for directory in which to place intermediate files
# Takes the number of lines of each of the two headers
# Adapted from Li Lei's Tutorial
#!/bin/env bash
#08/03/2018
# Fixes flip and force errors so that two vcf files which have already been sorted and filtered so that they have the same SNPs 
# will also have the same major and minor allele assignments and no missing values.
# Takes two vcf files.
# Adapted from Li Lei's Tutorial
#set -e
#set -u
#set -o pipefail
module load plink/1.90b
module load vcftools_ML/0.1.14
vcf1=$1
vcf2=$2
#run vcf-tools to find concordances and discordances. Puts them in diff_test_sortedwild_9k.diff.sites_in_files
vcftools --vcf $1 --diff $2  --diff-site --out diff_test
#Step 3 in tutorial
#Find the SNPs that are discordant and put a list of them in 
awk '$4=="O"{print $0}' diff_test.diff.sites_in_files >diff_test
#Step4 in tutorial
#Get the flipped list
grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$6] && $7 == complement[$8]{print $0}' diff_test|cut -f1,2) $2|cut -f 3 >flipped_SNP_list
#First round: flip the strand
plink --vcf $1 --allow-extra-chr --flip flipped_SNP_list --keep-allele-order --recode vcf --out flipped
#Compare the flipped vcf file with 9k vcf file
vcftools --vcf flipped.vcf --diff $2 --diff-site --out flipped
#Step 5 in tutorial
#Get the forced reference list
grep -f <(awk -F"\t" -v OFS="\t" '$5 == $8 && $6 == $7 { print $0 }' <(grep "O" flipped.diff.sites_in_files)|cut -f1,2) $2|cut -f3,4 >forced_ref_alleles
#Force reference
plink --vcf flipped.vcf --allow-extra-chr --a2-allele forced_ref_alleles --keep-allele-order --recode vcf --out forced_flipped
#Compare forced reference to 9k file
vcftools --vcf forced_flipped.vcf --diff $2 --diff-site --out forced_flipped
awk '$4=="O"{print $0}' forced_flipped.diff.sites_in_files >forced_flipped
#Step 6 in tutorial
#Get the flipped list
grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$8] && $6 == complement[$7]{print $0}' <(awk '$4=="O"{print $0}' forced_flipped.diff.sites_in_files)|cut -f1,2) $2|cut -f3 >flipped_SNP_for_forced_ref_flipped
#Get the forced reference list.
grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$8] && $6 == complement[$7]{print $0}' <(awk '$4=="O"{print $0}' forced_flipped.diff.sites_in_files)|cut -f1,2) $2|cut -f3,4 >forced_ref_for_forced_ref_flipped
#Flipped first and then force reference
plink --vcf forced_flipped.vcf --allow-extra-chr --flip flipped_SNP_for_forced_ref_flipped --a2-allele forced_ref_for_forced_ref_flipped --keep-allele-order --recode vcf --out forced_flipped_forced_ref_flipped
#Compare the forced reference to the original vcf file
vcftools --vcf forced_flipped_forced_ref_flipped.vcf --diff $2 --diff-site --out forced_flipped_forced_ref_flipped
#End of tutorial
#Remove those SNPs with a missing allele still left over after fixing flips and forces by first sorting both and then running grep
#Gets position of missing snps
awk '$4=="O"{print $0}' forced_flipped_forced_ref_flipped.diff.sites_in_files>fffrf1
awk '{print $1, $2}' fffrf1>fffrf
#Sorts forced_flipped_forced_ref_flipped.vcf
vcf-sort forced_flipped_forced_ref_flipped.vcf >fffrf.vcf
#Removes missing SNPs
grep -vf fffrf fffrf.vcf >filteredfffrf.vcf

