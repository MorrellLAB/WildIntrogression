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
                key_name = tmp[0] + '_' + tmp[1]    # ID1 and ID2 columns
                rest = [tmp[2], tmp[3], tmp[4], tmp[5], tmp[6]]
                if key_name in introg_dict.keys():
                    introg_dict[key_name].append(rest)
                else:
                    introg_dict[key_name] = [rest]
    return introg_dict


def main(file_list, threshold):
    d = {}
    read_file(input_file=file_list, introg_dict=d)

    tmp_d = {}
    for key in d.keys():
        if len(d[key]) == 1 and int(d[key][0][4]) < int(threshold):
            tmp_d[key] = d[key][0]

    for key in tmp_d.keys():
        d.pop(key, None)

    for key in d.keys():
        for i in range(0, len(d[key])):
            print('\t'.join(key.split('_')), *d[key][i], sep='\t')


main(sys.argv[1], sys.argv[2])

