#!/bin/env bash
# Takes map, vcf and hmp.txt file, trims them of unshared positions, converts them to a ped file, after removing the first position, which for some reason prevents my poor program from running then runs various plink analysis tools.
set -e
set -o pipefail
set -u
module load plink/1.90b
module load python3
map=$1
vcf=$2
hmp=$3
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
# removes first line
head -n 1 common3.list > remove.txt
sort remove.txt > sortedremove.txt
comm -23 common3.list sortedremove.txt > common4.list
# Creates new lists with files with extra positions removed
grep -Fwf common4.list $1 > finished.map
grep -Fwf common4.list $2 > almost.vcf
grep -Fwf common4.list $3 > almost.hmp.txt
# adds back headers
cat vcfheader.txt almost.vcf > finished.vcf
cat hmpheader.txt almost.hmp.txt > finished.hmp.txt
# Runs a program which makes a ped file
python3 /panfs/roc/groups/9/morrellp/depie014/hmp_to_PLINK.py finished.hmp.txt finished.map finished.vcf > finished.ped
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
rm common4.list
rm almost.vcf
rm almost.hmp.txt
rm remove.txt
rm sortedremove.txt
#Creates bed, bim  and fam files with plink
plink --ped finished.ped --map finished.map --allow-extra-chr
#renames output
mv plink.bed finished.bed
mv plink.bim finished.bim
mv plink.fam finished.fam
# runs plink commands to generate useful statistics
plink --bfile finished  --freq --allow-extra-chr
plink --bfile finished --distance square --allow-extra-chr
plink --bfile finished --distance square ibs --allow-extra-chr
# Renames plink output files and puts them in a directory named Analysis. It moves the bed, bim, and fam files to a directory named Data
mv plink.dist finished.dist
mv plink.dist.id finished.dist.id
mv plink.frq finished.frq
mv plink.log finished.log
mv plink.mibs finished.mibs
mv plink.mibs.id finished.mibs.id
mv plink.nosex finished.nosex
mkdir Data
mv finished.map ./Data
mv finished.vcf ./Data
mv finished.hmp.txt ./Data
mv finished.bed ./Data
mv finished.bim ./Data
mv finished.fam ./Data
mv finished.ped ./Data
mkdir Analysis
mv finished.* ./Analysis
