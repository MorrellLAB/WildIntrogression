#!/bin/env bash
#Connor Depies
#10/30/2017
#Program to split a vcf file into 100 snp intervals sliding by 25 snps and analyze their IBD with plink
#1) sorted VCF file
# Read the file in parameter and fill the array named "array"
# from https://stackoverflow.com/questions/20294918/extract-file-contents-into-array-using-bash
getArray() {
    array=() # Create array
    while IFS= read -r line # Read a line
    do
        array+=("$line") # Append line to the array
    done < "$1"
}

module load plink/1.90b
vcf=$1
# Splits the VCF file by chromosome, taking out the chromosome number and SNP number from columns 1 and 3, and removing header lines
awk '{print $1"\t" $3}' $1 |grep -v "#">/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt
grep "chr1H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr1.txt
grep "chr2H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr2.txt
grep "chr3H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr3.txt
grep "chr4H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr4.txt
grep "chr5H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr5.txt
grep "chr6H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr6.txt
grep "chr7H" /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr7.txt
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr1.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr1.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr2.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr2.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr3.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr3.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr4.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr4.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr5.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr5.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr6.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr6.list
awk '{print $2}' /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr7.txt >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/chr7.list
# Make an array of the SNPs in each Chromosome
chrarray=(chr1.list chr2.list chr3.list chr4.list chr5.list chr6.list chr7.list)
# Iterate over this array
for i in ${chrarray[@]};
do
    # Create an array containing SNP names from each Chromosome
    getArray /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/$i
    # Find IBD of each Chromosome
    # Makes a frqx file which contains allele frequencies. We then use this in the --genome command to correct for the fact that our samples are not in Hardy-Weinberg equilibrium.
    plink --vcf $1 --allow-extra-chr --extract /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/${i} --freqx --out /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}_out
    plink --vcf $1 --extract /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/${i} --genome --geno .15 --read-freq /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}_out.frqx --out /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}_out  --allow-extra-chr;
    # Find length of the array
    len=${#array[@]}
    # Create sliding windows each containing 100 SNPS from a chromosome, and moving by 25 SNPs each iteration
    for j in `seq 1 25 $len`;
    do
        if (("$j+100"<"$len"))
        then 
            # source: https://stackoverflow.com/questions/20243467/write-bash-array-to-file-with-newlines
            # Write each of these lists of 100 SNP windows to a file
            
            echo ${array[@]:$j:100} >/panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/${i}${j}_extract.txt
            # Makes a frqx file which contains allele frequencies. We then use this in the --genome command to correct for the fact that our samples are not in HArdy-Weinberg equilibrium.
            plink --vcf $1 -extract /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/${i}${j}_extract.txt --freqx --out /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}${j}_out --allow-extra-chr
            # Find IBD of this interval
            plink --vcf $1 --extract /panfs/roc/groups/9/morrellp/depie014/introgression/SNPlistsfor318wild/${i}${j}_extract.txt --genome --geno .15 --read-freq /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}${j}_out.frqx --out /panfs/roc/groups/9/morrellp/depie014/introgression/IBD2/fullgenomefile/${i}${j}_out --allow-extra-chr
        fi 
    done
done

            


    
