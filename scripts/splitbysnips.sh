#!/bin/env bash
#Connor Depies
#10/30/2017
#Program to split a vcf file into 25 snp intervals and analyze their IBD with plink
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
# Split the VCF file by chromosome
awk '{print $1, $3}' $1 >chr.txt
grep "chr1H" chr.txt >chr1.txt
grep "chr2H" chr.txt >chr2.txt
grep "chr3H" chr.txt >chr3.txt
grep "chr4H" chr.txt >chr4.txt
grep "chr5H" chr.txt >chr5.txt
grep "chr6H" chr.txt >chr6.txt
grep "chr7H" chr.txt >chr7.txt
awk '{print $2}' chr1.txt >chr1.list
awk '{print $2}' chr2.txt >chr2.list
awk '{print $2}' chr3.txt >chr3.list
awk '{print $2}' chr4.txt >chr4.list
awk '{print $2}' chr5.txt >chr5.list
awk '{print $2}' chr6.txt >chr6.list
awk '{print $2}' chr7.txt >chr7.list
# Make an array of the SNPs in each Chromosome
chrarray=(chr1.list chr2.list chr3.list chr4.list chr5.list chr6.list chr7.list)
# Iterate over this array
for i in ${chrarray[@]};
do
    # Create an array containing SNP names from each Chromosome
    getArray $i
     # Find IBD of each Chromosome
    plink --vcf $1 --extract ${i} --genome --geno .15 --out ${i}_out  --allow-extra-chr;
    # Find length of the array
    len=${#array[@]}
    # Create sliding windows each containing 100 SNPS from a chromosome, and moving by 25 SNPs each iteration
    for j in `seq 1 25 $len`;
    do
        if (("$j+100"<"$len"))
        then 
            # source: https://stackoverflow.com/questions/20243467/write-bash-array-to-file-with-newlines
            # Write each of these lists of 100 SNP windows to a file
            echo ${array[@]:$j:100} >${i}${j}_extract.txt
            # Find IBD of this interval
            plink --vcf $1 --extract ${i}${j}_extract.txt --genome --geno .15 --out ${i}${j}_out --allow-extra-chr
        fi 
    done
done

            


    