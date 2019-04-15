#!/usr/bin/env python3

#   Chaochih Liu - May 10, 2018

import sys

def read_file(input_file, introg_dict):
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('Ind1'):
                continue
            else:
                tmp = line.strip().split('\t')
                key_name = tmp[0] + '_' + tmp[1] + '_' + tmp[2] # ID1, ID2, and Chr columns
                rest = [tmp[3], tmp[4], tmp[5], tmp[6]]
                if key_name in introg_dict.keys():
                    introg_dict[key_name].append(rest)
                else:
                    introg_dict[key_name] = [rest]
    return introg_dict


def main(file_list, threshold):
    d = {} # Initialize empty dictionary
    read_file(input_file=file_list, introg_dict=d)
    #   combining overlapping intervals, so we can calculate
    #   size of continuous interval
    #   necessary b/c for same chr, some overlap and some don't
    merged_d = {}
    for key in d.keys():
        base = d[key][0]
        base_index = 0
        # Add base key value pair to merged_d
        merged_win = [
            base[0].split('-')[0], # win_start
            base[0].split('-')[1], # win_end
            base[1], # start_pos
            base[2], # end_pos
        ]
        merged_d[key] = [merged_win]
        for i in range(0, len(d[key])):
            if len(d[key]) == 1:
                merged_win = [
                    base[0].split('-')[0], # win_start
                    d[key][i][0].split('-')[1], # win_end
                    base[1], # start_pos
                    d[key][i][2], # end_pos
                ]
                merged_d[key] = [merged_win]
            elif (i > base_index and int(d[key][i][0].split('-')[0]) < int(d[key][i-1][0].split('-')[1])):
                # Values will contain: win_start, win_end, start_pos, end_pos
                merged_win = [
                    base[0].split('-')[0], # win_start
                    d[key][i][0].split('-')[1], # win_end
                    base[1], # start_pos
                    d[key][i][2], # end_pos
                ]
                # Dictionary values are a list of lists
                # Get last index of list of lists
                # This works because, other condition appends a list to list of lists
                # So, most current one we are updating will always be the nth list
                # If dictionary is not empty, add by list index
                if (key in merged_d.keys() and len(merged_d[key]) == 1):
                    merged_d[key][0] = merged_win
                elif (key in merged_d.keys() and len(merged_d[key]) > 1):
                    last_list = len(merged_d[key]) - 1
                    merged_d[key][last_list] = merged_win
                else:
                    merged_d[key] = [merged_win]
            else:
                if int(d[key][i][0].split('-')[0]) > int(d[key][i-1][0].split('-')[1]):
                    # Update base case
                    base = d[key][i]
                    base_index = i
                    # Values will contain: win_start, win_end, start_pos, end_pos
                    merged_win = [
                        base[0].split('-')[0], # win_start
                        d[key][i][0].split('-')[1], # win_end
                        base[1], # start_pos
                        d[key][i][2], # end_pos
                    ]
                    merged_d[key].append(merged_win)

    #   Calculate size of window
    for key in merged_d.keys():
        for i in range(0, len(merged_d[key])):
            #   Get total physical size of window
            phys_size = int(merged_d[key][i][3]) - int(merged_d[key][i][2])
            merged_d[key][i].append(str(phys_size))

    #   Print header line
    print("Ind1\tInd2\tChr\tSNP_Win_Start\tSNP_Win_End\tPhysPos_Start\tPhysPos_End\tInt_Phys_Size")
    #   Filter out interval sizes smaller than threshold
    #   Keep only those above threshold
    for key in merged_d.keys():
        for i in range(0, len(merged_d[key])):
            if int(merged_d[key][i][4]) >= int(threshold):
                print('\t'.join(key.split('_')), *merged_d[key][i], sep='\t')


main(sys.argv[1], sys.argv[2])
