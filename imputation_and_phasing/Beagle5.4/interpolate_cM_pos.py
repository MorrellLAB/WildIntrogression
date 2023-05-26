#!/usr/bin/env python3
"""Interpolate the genetic map positions of SNPs with unknown positions
in the PLINK 1.9 MAP file. cM positions will be interpolated based on
physical position order. SNPs where cM order disagree with physical position
order will be corrected based on the physical positions of flanking SNPs. This
assumes no physical positions are missing.

Usage: ./interpolate_cM_pos.py [plink.sorted_by_phys.map]

Where:
1) [plink.sorted_by_phys.map] is a Plink 1.9 map file sorted by the physical positions column.
Example: sort -k1,1 -k4,4n plink.map > plink.sorted_by_phys.map

Problems to resolve:
1) cM order disagrees with physical position order - correct cM position
    Need to correct order first, otherwise interpoloation for missing cM is incorrect
2) missing cM position in genetic map - interpolate position
3) Duplicate markers where cM AND physical positions are the same for different marker IDs - offset cM by 0.001
"""

import sys
import os
import pandas as pd

def parse_map(mapfile):
    """Return the MAP data as a dict with chromosome as the key."""
    map_data = []
    with open(mapfile, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            map_data.append([tmp[0], tmp[1], float(tmp[2]), int(tmp[3])])
            # if tmp[0] in map_dict.keys():
            #     map_dict[tmp[0]].append(tmp)
            # else:
            #     map_dict[tmp[0]] = [tmp[0], tmp[1], float(tmp[2], int(tmp[3]))]
    return map_data


map_fp = os.path.expanduser("~/Downloads/temp_msi/wbdc_bopa_snps.polymorphic.filt_miss_het.sorted_by_phys.map")

# Load map file
map_data = parse_map(map_fp)
# Convert to pandas data frame
map_df = pd.DataFrame(map_data, columns=['chr', 'id', 'cM', 'phys_pos'])
# Pull SNPs with non-missing cM locations
map_nonzero_cM = map_df[map_df['cM'] > 0].sort_values(by=['chr', 'phys_pos']).values.tolist()
# Convert to dict for easier processing by chromosome
map_nonzero_cM_dict = {}
for elem in map_nonzero_cM:
    if elem[0] in map_nonzero_cM_dict.keys():
        map_nonzero_cM_dict[elem[0]].append(elem)
    else:
        map_nonzero_cM_dict[elem[0]] = [elem[0:4]]

# Check if cM order agrees with physical position order
# If cM order disagrees with physical order, interpolate cM location based on physical positions
map_corrected_order = []
for chr_key in map_nonzero_cM_dict.keys():
    for idx, elem in enumerate(map_nonzero_cM_dict[chr_key]):
        last_idx = len(map_nonzero_cM_dict[chr_key]) - 1
        if idx == 0:
            downstream_cM = map_nonzero_cM_dict[chr_key][idx+1][2]
            if elem[2] <= downstream_cM:
                # Order is correct, nothing needs to be done
                map_corrected_order.append(elem)
            else:
                # Use closest downstream variant's cM position since there's no midpoint that can be calculated
                corrected_cM = downstream_cM
                map_corrected_order.append([elem[0], elem[1], corrected_cM, elem[3]])
        elif idx != 0 and idx != last_idx:
            upstream_cM = map_nonzero_cM_dict[chr_key][idx-1][2]
            downstream_cM = map_nonzero_cM_dict[chr_key][idx+1][2]
            if upstream_cM <= elem[2] <= downstream_cM:
                # Order is correct, nothing neesd to be done
                map_corrected_order.append(elem)
            else:
                # Correct cM position to be midpoint between upstream/downstream marker
                # But, we could have two consecutive downstream markers where both cM
                # locations have incorrect order, so find next SNP that has cM >= current cM position
                print(idx, elem)
                if elem[2] > downstream_cM:
                    increment_idx = 1
                    #!!!!!!! This while loop isn't working as expected to keep searching for the right SNP to use to interpolate CM position
                    while downstream_cM <= elem[2]:
                        # Check if we are at the last index in our list
                        if idx+increment_idx <= last_idx:
                            print(idx + increment_idx)
                            downstream_cM = map_nonzero_cM_dict[chr_key][idx+increment_idx][2]
                            increment_idx += 1
                        elif idx+increment_idx > last_idx:
                            # If we are at the last index in our list and the cM is still less than our current elem, leave cM position as is and "reset" downstream_cM
                            map_corrected_order.append(elem)
                            downstream_cM = map_nonzero_cM_dict[chr_key][idx+1][2]
                if upstream_cM > elem[2]:
                    decrement_idx = 1
                    # Note: this is the updated downstream_cM
                    interpolated_cM = (upstream_cM + downstream_cM) / 2
                    while interpolated_cM < map_nonzero_cM_dict[chr_key][idx-1-decrement_idx][2]:
                        # Re-calculate interpolated
                        upstream_cM = map_nonzero_cM_dict[chr_key][idx-1-decrement_idx][2]
                        interpolated_cM = (upstream_cM + downstream_cM) / 2
                # Final calculation of interpolated_cM
                interpolated_cM = (upstream_cM + downstream_cM) / 2
                map_corrected_order.append([elem[0], elem[1], interpolated_cM, elem[3]])
        elif idx == last_idx:
            upstream_cM = map_nonzero_cM_dict[chr_key][idx-1][2]
            if upstream_cM <= elem[2]:
                # Order is correct, nothing needs to be done
                map_corrected_order.append(elem)
            else:
                # Use closest upstream variant's cM position since there's no midpoint that can be calculated
                corrected_cM = upstream_cM
                map_corrected_order.append([elem[0], elem[1], corrected_cM, elem[3]])


# Store interpolated cM positions
new_cM_map = []
for chr_key in map_data.keys():
    last_idx = len(map_data[chr_key]) - 1
    for idx, elem in enumerate(map_data[chr_key]):
        if idx == 0:
            # No upstream marker available for interpolation
            new_cM_map.append(elem)
        elif idx == last_idx and elem[2] == '0':
            # Last index in list for current chrom, no more downstream SNPs
            # and cM position is missing
            # Use upstream SNP's cM position since we can't get a midpoint
            cM_interpolated = map_data[chr_key][idx-1][2]
            new_cM_map.append([elem[0], elem[1], str(cM_interpolated), elem[3]])
        elif idx != 0 and idx != last_idx:
            if elem[2] == '0' and idx != last_idx:
                # Starting increment/decrement
                decrement_idx = 1
                increment_idx = 1
                # interpolate cM location based on physical position order
                upstream_elem = map_data[chr_key][idx-decrement_idx]
                downstream_elem = map_data[chr_key][idx+increment_idx]
                # Continue until we find closest SNP with non-zero cM
                while upstream_elem[2] == '0':
                    decrement_idx += 1
                    upstream_elem = map_data[chr_key][idx-decrement_idx]
                while downstream_elem[2] == '0':
                    increment_idx += 1
                    downstream_elem = map_data[chr_key][idx+increment_idx]
                # Get new cM position using midpoint between SNPs with known cM
                cM_interpolated = (float(upstream_elem[2]) + float(downstream_elem[2])) / 2
                # Add interpolated cM
                new_cM_map.append([elem[0], elem[1], str(cM_interpolated), elem[3]])
            else:
                 new_cM_map.append(elem)
