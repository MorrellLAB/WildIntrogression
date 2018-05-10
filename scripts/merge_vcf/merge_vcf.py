#!/usr/bin/env python3

#   Chaochih Liu - April 16, 2018

"""This script takes in two VCF sorted (numerically by physical position) files and merges them.
Important Caveat:
    Output VCF file is unsorted and NOT de-duplicated. Sorting within this script and de-duplicating
    feature will be added in the future. For now, we will have to run the following extra Linux command line after:

    (head -n 18 test_merged.vcf && tail -n +19 test_merged.vcf | sort -uV -k1,2) > test_merged_sorted.vcf

Usage: ./merge_vcf.py [vcf_file1.vcf] [vcf_file2.vcf] [contig_file] > out_file.vcf

Where:
    [file1.vcf] and [file2.vcf] are sorted by CHROM and then POS columns.
    [contig_file] is a file that contains correct contig length in the '##contig' format
"""

import sys

def read_meta_info(vcf_file):
    """Function that reads in VCF file meta-information ('##') lines only and stores each line in a list."""
    meta_info_dat = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                meta_info_dat.append(line.strip())
            else:
                break
    return meta_info_dat


def read_vcf(vcf_file):
    """Function that reads in VCF file, excluding meta-information lines and including header line starting with #CHROM.
    This stores each row in a dictionary in the following format:
    '#CHROM_POS': [[CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT], [CIho10033, CIho10034, ...]]"""
    vcf_dict = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            else:
                tmp = line.strip().split('\t')
                key_name = tmp[0] + '_' + tmp[1] # 0 - #CHROM; 1 - POS
                #   Mandatory columns according to VCF 4.1 specification
                #   except FORMAT column is only present when genotype data is present
                mandatory_cols = tmp[0:9] # Pulls out columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
                geno = tmp[9:] # Pulls out all columns with genotype data (i.e. '0/0', '1/1', etc)
                vcf_dict[key_name] = [mandatory_cols, geno]
    return vcf_dict


def merge_meta_info(contig_pos, meta1, meta2):
    """Function that merges VCF meta-info fields and replaces Plink's contig positions with WBDC exon capture contig positions.

    Plink uses the last SNP in each chromosome to calculate the length of each chromosome for the contig meta-information field.
    The length of the chromosomes will change depending on the SNPs present in the VCF file, so it is not accurate.
    We want to replace the ##contig meta-info field with the actual length of each chromosome, which can be found
    in the WBDC exon capture VCF file.

    Required meta-info line for standard VCF file: ##fileformat=VCFv4.2

    IMPORTANT: this function currently only handles the following meta-info lines:
    - ##fileformat
    - ##FORMAT
    - ##contig

    To handle additional meta-info lines, code needs to be modified!
    """
    c = meta1 + meta2   # concatenate lists
    u = set(c)  # use set to pull out unique elements in list
    m = []  # initialize empty list
    for line in u:
        if '##fileformat' in line:
            m.insert(0, line)
        elif '##FORMAT' in line:
            m.append(line)

    #   Now add correct contig lengths from exon capture VCF files
    for line in contig_pos:
        m.append(line)
    return m


def merge_vcf(d1, d2):
    """Function that merges two VCF files based on matching dictionary keys. Takes in two VCF files stored in dictionaries"""
    for key in d1.keys():
        if key in d2.keys():
            print('\t'.join(d2[key][0] + d1[key][1] + d2[key][1]))
        elif key not in d2.keys():
            print('\t'.join(d1[key][0] + d1[key][1] + (['./.'] * len(d2['#CHROM_POS'][1]))))
            for key2 in d2.keys():
                if key2 not in d1.keys():
                    print('\t'.join(d2[key2][0] + (['./.'] * len(d1['#CHROM_POS'][1])) + d2[key2][1]))


def main(vcf_file1, vcf_file2, contig_file):
    """Function that runs the program."""
    #   Read in meta-info lines
    minfo1=read_meta_info(vcf_file=vcf_file1)
    minfo2=read_meta_info(vcf_file=vcf_file2)
    contig=read_meta_info(vcf_file=contig_file)

    #   Read in vcf files
    v1=read_vcf(vcf_file=vcf_file1)
    v2=read_vcf(vcf_file=vcf_file2)

    #   Merge meta-info lines
    merged_header=merge_meta_info(contig_pos=contig, meta1=minfo1, meta2=minfo2)
    for line in merged_header:
        print(line)

    #   Merge remaining vcf file
    merge_vcf(d1=v1, d2=v2)


main(sys.argv[1], sys.argv[2], sys.argv[3]) # Run the program
