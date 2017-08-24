#!/bin/env bash
#Connor Depies Aug 17, 2017
#This program takes 4 inputs:
#1 ) A map file
#2 ) A partial vcf file
#3 ) A genotype file
#These three files are filtered of unshared snps, then converted to a full vcf file. Then the program fixes flip and force errors, removes snps which don't have both alleles.
#4 ) A full vcf file
#We then merge are full vcf files
#Connor Depies August, 18, 2017
set -e
set -o pipefail
set -u
module load plink/1.90b
module load python3
module load vcftools_ML/0.1.14
map=$1
vcf=$2
hmp=$3
vcf2=$4
# Stores headers
head -n 8 $2  > vcfheader.txt 
# Make sure you have 8 lines in your vcfheader
head -n 1 $3 > hmpheader.txt
# Removes positions not in all three input files
awk '{print $2}' $1 > map.list
awk '{print $3}' $2 > vcf.list
awk '{print $1}' $3 > hmp.list
sort map.list > sortedmap.list
sort vcf.list >sortedvcf.list
sort hmp.list >sortedhmp.list
comm -12 sortedvcf.list sortedmap.list > common1.list
comm -12 sortedvcf.list sortedhmp.list > common2.list
comm -12 common1.list common2.list > common3.list
# Creates new lists with files with extra positions removed
grep -Fwf common3.list $1 > diogenes.map
grep -Fwf common3.list $2 > almost.vcf
grep -Fwf common3.list $3 > almost.hmp.txt
# adds back headers
cat vcfheader.txt almost.vcf > diogenes.vcf
cat hmpheader.txt almost.hmp.txt > diogenes.hmp.txt
# Runs a program which makes a ped file
python3 /panfs/roc/groups/9/morrellp/depie014/hmp_to_PLINK.py diogenes.hmp.txt diogenes.map diogenes.vcf > diogenes.ped
# removes unnecessary files
rm vcfheader.txt
rm hmpheader.txt
rm map.list
rm vcf.list
rm hmp.list
rm sortedmap.list
rm sortedvcf.list
rm sortedhmp.list
rm common1.list
rm common2.list
rm common3.list
rm almost.vcf
rm almost.hmp.txt
#Creates bed, bim  and fam files with plink
plink --ped diogenes.ped --map diogenes.map --allow-extra-chr
#renames output
mv plink.bed diogenes.bed
mv plink.bim diogenes.bim
mv plink.fam diogenes.fam
#Makes new vcf_file
plink --bfile diogenes --recode vcf-iid --out wild_9k --allow-extra-chr
#Fixes flips and forces
#Sort vcf file by position
#step 2 in tutorial
vcf-sort wild_9k.vcf>sortedwild_9k.vcf
#run vcf-tools to find concordances and discordances. Puts them in diff_test_sortedwild_9k.diff.sites_in_files
vcftools --vcf sortedwild_9k.vcf --diff $2  --diff-site --out diff_test
#Step 3 in tutorial
#Find the SNPs that are discordant and put a list of them in 
awk '$4=="O"{print $0}' diff_test.diff.sites_in_files >diff_test
#Step4 in tutorial
#Get the flipped list
grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$6] && $7 == complement[$8]{print $0}' diff_test|cut -f1,2) $2|cut -f 3 >flipped_SNP_list
#First round: flip the strand
plink --vcf wild_9k.vcf --allow-extra-chr --flip flipped_SNP_list --keep-allele-order --recode vcf --out flipped
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
#Step 6 yields empty files, since there are no forces left to fix, and an empty file is causing a problem in running the program, so I will deal
#Step 6 in tutorial
#Get the flipped list
#grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$8] && $6 == complement[$7]{print $0}' <(awk '$4=="O"{print $0}' forced_ref_flipped.diff.sites_in_files)|cut -f1,2) $2|cut -f3 >flipped_SNP_for_forced_ref_flipped
#Get the forced reference list.
#grep -f <(awk -F"\t" -v OFS="\t" 'BEGIN {complement["A"] = "T"; complement["T"] = "A"; complement["C"] = "G"; complement["G"] = "C";} $5 == complement[$8] && $6 == complement[$7]{print $0}' <(awk '$4=="O"{print $0}' forced_ref_flipped.diff.sites_in_files)|cut -f1,2) $2|cut -f3,4 >forced_ref_for_forced_ref_flipped
#Flipped first and then force reference
#plink --vcf forced_ref_flipped.vcf --allow-extra-chr --flip flipped_SNP_for_forced_ref_flipped --a2-allele forced_ref_for_forced_ref_flipped --keep-allele-order --recode vcf --out forced_flipped_forced_ref_flipped
#Compare the forced reference to the original vcf file
#vcftools --vcf forced_flipped_forced_ref_flipped.vcf --diff $2 --diff-site --out forced_flipped_forced_ref_flipped
#Remove those SNPs with a missing allele still left over after fixing flips and forces by first sorting both and then running grep
vcf-sort forced_flipped >sorted_forced_flipped
#Gets position of missing snps
awk '{print $2}' forced_flipped>ff
vcf-sort forced_flipped.vcf >sorted_forced_flipped.vcf
#Removes them
grep -vf ff sorted_forced_flipped.vcf >filteredsff.vcf
#Filter out SNPs not shared by both filteredsff.vcf and $4
#Removes and stores headers (filteredsff header is 13 lines long, $4 is also 13 lines long)
head -n 13 filteredsff.vcf >fsffheader.txt
tail -n +14 filteredsff.vcf >h_fsff.vcf
head -n 13 $4 >NAMheader.txt
tail -n +14 $4 >h_NAM.vcf
# Makes uncommon list.
awk '{print $3}' h_fsff.vcf >h_fsff.list
awk '{print $3}' h_NAM.vcf >h_NAM.list
sort h_fsff.list >sh_fsff.list
sort h_NAM.list >sh_NAM.list
comm -12 sh_fsff.list sh_NAM.list >fstt_NAM.list
#Filters
grep -Fwf fstt_NAM.list h_fsff.vcf >filtered_NAM.vcf
grep -Fwf fstt_NAM.list h_NAM.vcf >filtered_fstt.vcf
#Filter out snps shared by both which are missing alleles.
#Get lists of those positions which are missing an allele.
grep -v "#" filtered_fstt.vcf  | cut -f 1,2,4,5>ffstt.list
grep -v "#" filtered_NAM.vcf  | cut -f 1,2,4,5>NAM.list
#Get list of all missing positions
diff -y ffstt.list NAM.list| grep '|' | cut -f 1,2 >missing_both
echo "Walrus"
#Do the filtering
grep -vf missing_both filtered_fstt.vcf >filteredmissingffstt.vcf
grep -vf missing_both filtered_NAM.vcf >filteredmissingNAM.vcf
#Restore Headers
cat fsffheader.txt filteredmissingffstt.vcf >hfilteredmissingffstt.vcf
cat NAMheader.txt filteredmissingNAM.vcf >hfilteredmissingNAM.vcf
#Sort
vcf-sort hfilteredmissingffstt.vcf >finalwild_9k.vcf
vcf-sort hfilteredmissingNAM.vcf >finalNAM.vcf
#Zip vcf files
bgzip -c finalwild_9k.vcf >finalwild_9k.vcf.gz
bgzip -c finalNAM.vcf >finalNAM.vcf.gz
