#!/usr/bin/env python3
"""Adds the cM positions to the Plink MAP file. If SNP exists in .map file, but
does not have a position in the consensus genetic map, exclude the SNP. If SNP ID in
MAP file doesn't have a cM position, the SNP ID will be saved in a separate output file
in the same directory as the MAP file with file prefix "no_cM_pos_id".

Usage: ./add_cM_to_MAP.py [map_file] [cM_pos_file] [cm_column_number] [variant_id_column_number] > updated.map

Where:
1) [map_file] is the filepath to the Plink 1.9 MAP file
2) [cM_pos_file] is the filepath to the a tab-delimited file containing
    the columns with SNP ID and cM position.
3) [cm_column_number] is the column number (1-based indexing) that has the cM position.
4) [variant_id_column_number] is the column number (1-based indexing) that has the variant ID.

Note: We will subtract 1 from the column numbers provided to get the list index
since Python uses 0-based indexing.
"""

import sys
import os

# Print usage message if no user arguments are provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def read_map(cm_pos_file, id_index):
    """Read in full pedigree CSV file and store in a dictionary.
    This file contains info on the Family ID, Maternal ID, and Paternal ID"""
    cm_pos = {}
    with open(cm_pos_file, 'r') as f:
        for idx, elem in enumerate(f):
            tmp = elem.strip().split('\t')
            snpid = tmp[id_index]
            # Add to dictionary
            cm_pos[snpid] = tmp
    return(cm_pos)


def main(map_file, cm_pos_file, cm_col_num, id_col_num):
    """Read through Plink PED file and add family information."""
    # Prepare array indices
    cm_index = int(cm_col_num) - 1
    id_index = int(id_col_num) - 1
    # Read in map file
    cm_dict = read_map(os.path.expanduser(cm_pos_file), id_index)
    # Read in MAP file and add the cM positions
    #   The third column is the position in morgans or centimorgans
    # Keep track of SNPs in MAP file that don't have cM position
    no_cm_pos = []
    with open(os.path.expanduser(map_file), 'r') as f:
        for line in f:
            tmp_map = line.strip().split()
            curr_snpid = tmp_map[1]
            if curr_snpid in cm_dict.keys():
                # Replace the missing cM position
                tmp_map[2] = cm_dict[curr_snpid][cm_index]
                # Print line
                print('\t'.join(tmp_map))
            else:
                no_cm_pos.append(curr_snpid)
    # Save no cM position variant IDs to separate file
    map_file_dir = os.path.dirname(os.path.abspath(map_file))
    map_file_bn = os.path.splitext(os.path.basename(os.path.abspath(map_file)))[0]
    out_fp = map_file_dir + '/no_cM_pos_ids-' + map_file_bn + '.txt'
    with open(out_fp, 'w') as outfile:
        outfile.write('\n'.join(item for item in no_cm_pos))


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
