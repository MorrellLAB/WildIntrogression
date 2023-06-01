#!/usr/bin/env python3
"""This script prepares the pops.txt file required for
the ABBA BABA test by Simon Martin:
https://github.com/simonhmartin/genomics_general.
This script prints to stdout a pops.txt file with 2 columns:
1) sample name, 2) population label

Usage: ./generate_pops_file.py [samp_info_fp] [vcf_samp_fp] > pops.txt

Where:
1) [samp_info_fp] is a file where column 1 contains the accession ID and
    column 5 contains the population label.
2) [vcf_samp_fp] is a file that contains 1 column of sample names as they
    appear in the VCF or .geno.gz file."""

import sys
import os

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def read_sample_info(samp_info_fp):
    """Read CSV file that contains group labels
    for each accession ID."""
    sample_dict = {}
    with open(samp_info_fp, 'r') as f:
        for line in f:
            tmp = line.strip().split(',')
            sample_dict[tmp[0]] = tmp
    return(sample_dict)


def main(samp_info_fp, vcf_samp_fp):
    """Driver function."""
    samp_dict = read_sample_info(os.path.expanduser(samp_info_fp))
    with open(os.path.expanduser(vcf_samp_fp), 'r') as f:
        for line in f:
            tmp_vcf_samp = line.strip()
            if tmp_vcf_samp in samp_dict.keys():
                print('\t'.join([tmp_vcf_samp, samp_dict[tmp_vcf_samp][4]]))


main(sys.argv[1], sys.argv[2])
