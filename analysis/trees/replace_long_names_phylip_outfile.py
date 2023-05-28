#!/usr/bin/env python3
"""Phylip cuts names short if they are longer than 10 characters, this
is problematic for figtree because the accession names are no longer unique.
Replace names >10 characters with shortname in a lookup table."""

import sys
import os

# User provided input arguments
phylip_fp = os.path.expanduser(sys.argv[1])
lookup_fp = os.path.expanduser(sys.argv[2])

def read_lookup(lookup_fp):
    lookup = {}
    with open(lookup_fp, 'rt') as file:
        for line in file:
            tmp = line.strip().split('\t')
            lookup[tmp[0]] = tmp[1]
    return lookup

lookup = read_lookup(lookup_fp)

with open(phylip_fp, 'rt') as file:
    for line in file:
        tmp = line.strip().split()
        if tmp[0] in lookup.keys():
            # Replace name with short name in lookup
            print('\t'.join([lookup[tmp[0]], tmp[1]]))
        else:
            print(line.strip())
