#!/usr/bin/env python3
"""Convert Plink .fam file to Pong population labels file.

Usage: ./fam_to_pong_pop.py [plink.fam] [sample_codes_fp] > ind2pop.txt

Where:
1) [plink.fam] where column 1 contains the accession codes that correspond to column 2 in [sample_codes_fp]

2) [sample_codes_fp] is tab delimited and formatted as pop_name, pop_code. Example:
wild    1
wild_introgressed    2
landrace    3
cultivar    4
breeding    5
genetic    6
uncertain    7

Important to keep the order the same as in the .fam file since those correspond to rows in the Admixture output files.
"""

import sys
import os

if len(sys.argv) < 2:
    print(__doc__)
    exit(0)

# User provided input arguments
plink_fam_fp = os.path.expanduser(sys.argv[1])
sample_codes_fp = os.path.expanduser(sys.argv[2])

def read_acc_codes(sample_codes_fp):
    acc_codes = {}
    with open(sample_codes_fp, 'rt') as file:
        for line in file:
            tmp = line.strip().split()
            acc_codes[tmp[1]] = tmp[0]
    return acc_codes

# Load accession codes corresponding to names
acc_codes = read_acc_codes(sample_codes_fp)

# Translate code to names from column 1 in .fam file
with open(plink_fam_fp, 'rt') as file:
    for line in file:
        tmp = line.strip().split()
        print(acc_codes[tmp[0]])
