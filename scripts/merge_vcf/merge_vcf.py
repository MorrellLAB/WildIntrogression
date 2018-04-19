#!/usr/bin/env python3

#   Chaochih Liu - April 16, 2018

"""This script takes in two VCF sorted (numerically by physical position) files and merges them.

Usage: ./merge_vcf.py [file1.vcf] [file2.vcf] > out_file.vcf

Where:
    [file1.vcf] and [file2.vcf] are sorted by CHROM and then POS columns.
"""

import sys

def read_meta_info(vcf_file):
    """Function that reads in VCF file meta-information ('##') lines only."""
    meta_info_dat = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                meta_info_dat.append(line.strip().split('\n'))
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



