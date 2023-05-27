#!/usr/bin/env python3
"""Convert Annovar .exonic_variant_function to BAD_Mutations .subs format.

Assumes nonsynonymous variants have already been extracted and the exonic
annotation format is c.C122T and NOT HGVS format (c.122C>T rather than c.C122T).

The script currently handles the format c.C122T. For more info, see help message:
"annotate_variation.pl --version" and look for description under the flag --hgvs.
More info on sequence variant nomenclature can be found here: https://varnomen.hgvs.org/bg-material/simple/

Usage: ./ANNOVAR_to_subs.py [annovar.exonic_variant_function] [output_long_subs.txt] [out_dir]

Where:
1) [annovar.exonic_variant_function] is the Annovar output file .exonic_variant_function subset for nonsynonymous variants as defined by Sequence Ontology
2) [output_long_subs.txt] is the path to the output long_subs.txt file
3) [out_dir] is the path to the output directory that will store the .subs files grouped by transcript ID

Outputs the following files:

Output File 1 format: long_subs.txt
1) Transcript identifier
2) protein position
3) Reference amino acid
4) Variant amino acid
5) Variant identifier

Output File 2 format: .subs
1) protein position
2) unique variant identifier
"""

import sys
import os
import re
import pandas as pd

# Print usage message if no input arguments are provided
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

# User provided input argument
annovar_exonic_var = os.path.abspath(sys.argv[1])
long_subs_fp = os.path.abspath(sys.argv[2])
out_dir = os.path.abspath(sys.argv[3])

# Load file and store relevant fields in list
def read_annovar_exonic_var(annovar_exonic_var):
    annovar_list = []
    with open(annovar_exonic_var, 'rt') as file:
        for line in file:
            tmp = line.strip().split('\t')
            annovar_list.append([tmp[2], tmp[11], tmp[12], tmp[14], tmp[15]])
    return annovar_list

# Prepare long_subs.txt file
# Store in dictionary so that we can make sure that the substitutions are
# grouped by transcript before writing to .subs files
def reformat_subs(annovar_list):
    long_subs_list = []
    subs_dict = {}
    for elem in annovar_list:
        var_id = '_'.join([elem[1], elem[2]]) + '_' + '/'.join([elem[3], elem[4]])
        tmp_info = elem[0].split(':')
        transcript_id = tmp_info[1]
        protein_info = tmp_info[4]
        # Check if we are pulling the protein level info (should start with "p.")
        if protein_info.startswith("p."):
            # Reformat refAA, variantAA, and protein position info
            # Split string by letters and numbers
            # Example string format: p.A20T
            split_protein_info = re.split('(\d+)', protein_info.split('.')[1])
            ref_aa = split_protein_info[0]
            var_aa = split_protein_info[2].strip(',')
            protein_pos = split_protein_info[1]
            # Add to long subs list
            long_subs_list.append([transcript_id, protein_pos, ref_aa, var_aa, var_id])
            # Add to subs dict grouped by transcript id
            if transcript_id in subs_dict.keys():
                subs_dict[transcript_id].append([protein_pos, var_id])
            else:
                subs_dict[transcript_id] = [[protein_pos, var_id]]
        else:
            sys.stderr.write('Error: protein position info not available, please check if protein level info exists in the annovar_output.exonic_variant_function file (should be in column 3 and start with p.)\n')
            exit(1)
    return (long_subs_list, subs_dict)

# Convert annovar to subs
annovar_list = read_annovar_exonic_var(annovar_exonic_var)
long_subs_list, subs_dict = reformat_subs(annovar_list)
# Save long_subs to file
pd.DataFrame(long_subs_list).to_csv(long_subs_fp, sep="\t", header=False, index=False)
# Save transcripts to .subs files
for transcript in subs_dict.keys():
    subs_fp = out_dir + '/' + transcript + '.subs'
    pd.DataFrame(subs_dict[transcript]).to_csv(subs_fp, sep="\t", header=False, index=False)
