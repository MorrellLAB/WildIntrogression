#!/usr/bin/env python3

#   Chaochih Liu - May 10, 2018

import sys
from decimal import Decimal

def read_list_of_files(input_list):
    with open(input_list, 'r') as f:
        lines = f.read().splitlines()
    return lines


def read_file(input_file, introg_dict):
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('IID1'):
                continue
            else:
                tmp = line.strip().split('\t')
                key_name = tmp[0] + '_' + tmp[1]    # ID1 and ID2 columns
                pi_hat = Decimal(tmp[7].strip(''))
                if key_name in introg_dict.keys():
                    introg_dict[key_name].append(pi_hat)
                else:
                    introg_dict[key_name] = [pi_hat]
    return introg_dict


def main(file_list):
    #   Read in list of files
    wc_list = read_list_of_files(input_list=file_list)

    d = {}
    for i in wc_list:
        read_file(input_file=i, introg_dict=d)

    out_dict = {k:sum(v) for k, v in d.items()}

    for key in out_dict.keys():
        print('\t'.join(key.split('_')), '\t', out_dict[key])


main(sys.argv[1])