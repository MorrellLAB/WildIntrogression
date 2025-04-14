#!/usr/bin/env python3
"""Generate .pop file compatible with ADMIXTURE's supervised analysis.

Usage: ./plink_fam_to_pop_file.py file.fam pop_info.csv > outfile.pop
"""

import sys
import os

# Print usage message if there aren't any argument provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def load_pop_info(pop_info_fp):
    pop_dict = {}
    with open(os.path.expanduser(pop_info_fp), 'r') as f:
        for line in f:
            tmp=line.strip().split(',')
            samp_id=tmp[0]
            samp_group=tmp[4]
            pop_dict[samp_id] = samp_group
    return pop_dict


def main(plink_fam_fp, pop_info_fp):
    """Main driver function."""
    # Read in file containing population grouping for individuals with known vs unknown ancestry
    pop_info_dict = load_pop_info(pop_info_fp)
    # Read in plink .fam file and generate .pop file
    with open(os.path.expanduser(plink_fam_fp), 'r') as f:
        for line in f:
            tmp=line.strip().split(' ')
            fam_samp_id=tmp[1]
            # If key in dict, use ancestry group label in dict, else assume unknown ancestry
            if fam_samp_id in pop_info_dict:
                print(pop_info_dict[fam_samp_id])
            else:
                print('-')


main(sys.argv[1], sys.argv[2])
