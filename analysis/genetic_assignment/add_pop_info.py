#!/usr/bin/env python3
"""This script adds the population information to a .txt
file containing sample names.

Usage: ./add_pop_info.py [sample_list_fp] [pop_codes_fp] [pop_info_xlsx_fp] [my_sheet_name] > sample_names_with_pop_info.txt

Where:
1) [sample_list_fp] is a single column .txt file with one sample name per line. Samples
                    should be in the same order as they are in the original VCF file
2) [pop_codes_fp] is a two column file that has the population name (matches pop name
                    in [pop_info_xlsx_fp] sheet) and the corresponding numerical code
3) [pop_info_xlsx_fp] is an xlsx spreadsheet where at minimum there are the two columns,
                    "Accession_ID" and "Accession_Type". There can be additional columns
                    but those two required columns must be named as mentioned.
4) [my_sheet_name] is the name of the sheet in the xlsx file of interest.

Note: In the population info xlsx file, the population grouping
must be under a column named 'Accession_Type' otherwise this script
won't work.

Accessions without any population group info will be prefaced with an "NA" in the output.
This way we can go back and find the missing population information easily.
"""

import sys
import os
import pandas as pd

# Print usage message if there aren't any argument provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def main(sample_list_fp, pop_codes_fp, pop_info_xlsx_fp, my_sheet_name):
    """Main driver function."""
    # Read in xlsx file
    # In our xlsx file, the header row index is 1
    # And we have combined all pop info into a single sheet
    pop_df = pd.read_excel(os.path.expanduser(pop_info_xlsx_fp),
        sheet_name=my_sheet_name,
        header=1)

    # Read in file that contains population codes
    pop_codes = {}
    with open(os.path.expanduser(pop_codes_fp), 'r') as f:
        for line in f:
            tmp = line.strip().split()
            pop_codes[tmp[0]] = tmp[1]

    # Read in sample list and add population code to 2nd column
    with open(os.path.expanduser(sample_list_fp), 'r') as f:
        for line in f:
            curr_samp = line.strip()
            if pop_df['Accession_ID'].eq(curr_samp).any():
                # Extract row corresponding to sample name, then pull the
                #   Accession_Type from the row
                acc_type_pd = pop_df.loc[pop_df['Accession_ID'] == curr_samp, 'Accession_Type']
                acc_type = pd.Series.to_string(acc_type_pd).lower().split()[1]
                # Pull pop code
                print(':'.join([pop_codes[acc_type], curr_samp]))
            else:
                print(':'.join(["NA", curr_samp]))


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
