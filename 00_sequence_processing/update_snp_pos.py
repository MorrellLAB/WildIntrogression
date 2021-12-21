#!/usr/bin/env python3
"""Script that converts updates VCF SNP positions to Morex v3 positions and outputs
UNSORTED VCF files: 1) updated position vcf, and 2) missing SNPs vcf (if any)

Usage: ./update_snp_pos.py [vcf_file] [new_snp_pos_vcf] [out_dir] [out_filename]

Where:
1) [vcf_file] is the full filepath to the VCF file where we want to update the positions
2) [new_snp_pos_vcf] is a VCF file that contains the new positions
3) [out_dir] is the full filepath to the output directory
4) [out_filename] is the output filename of our vcf with updated positions
"""

import sys
import os

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)

vcf_file = "/Users/chaochih/Downloads/toy_wbdc_bopa.vcf"
new_snp_pos_vcf = "~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/bopa_idt95_noRescuedSNPs.vcf"


def read_vcf(vcf):
    """Read VCF file and store in a dictionary."""
    header_lines = []
    snps_dict = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                l = line.strip().split('\t')
                # Store in dictionary
                snps_dict[l[2]] = l
    return (header_lines, snps_dict)


def update_pos(vcf, new_pos_vcf):
    """Update SNP positions."""
    updated_dict = {}
    missing_snps = {}
    for key in vcf.keys():
        if key in new_pos_vcf.keys():
            # Key exists in new pos lookup vcf
            # New chromosome and positions
            new_chrom = new_pos_vcf[key][0]
            new_pos = new_pos_vcf[key][1]
            tmp_line = vcf[key]
            tmp_line[0] = new_chrom
            tmp_line[1] = new_pos
            updated_dict[key] = tmp_line
        else:
            # Key doesn't exist in new pos lookup vcf
            # Save SNP name to separate list
            missing_snps[key] = vcf[key]
    return (updated_dict, missing_snps)


def main(vcf_file, new_snp_pos_vcf, out_dir, out_filename):
    """Driver function."""
    # Read in vcf files
    vcf_header, vcf_dict = read_vcf(os.path.expanduser(vcf_file))
    new_pos_vcf_header, new_pos_vcf_dict = read_vcf(os.path.expanduser(new_snp_pos_vcf))

    # Update positions
    updated_vcf, missing_snp_names = update_pos(vcf_dict, new_pos_vcf_dict)

    # Save to files
    # Prepare output filepath
    out_fp = os.path.expanduser(out_dir.rstrip('/')) + "/" + out_filename
    # Start with a clean slate since we are appending
    if os.path.exists(out_fp):
        os.remove(out_fp)
    # Create and open file for writing
    # New VCF with update positions
    cof_new_vcf = open(out_fp, 'a+')
    for i in vcf_header:
        cof_new_vcf.write(i)
    for key in updated_vcf.keys():
        cof_new_vcf.write('\t'.join(updated_vcf[key]) + '\n')
    cof_new_vcf.close()
    # Check if we have missing SNPs
    # If so, save list of missing SNPs to separate VCF
    # "missing" means the SNP in vcf_file doesn't exist in new_snp_pos_vcf
    if len(missing_snp_names) > 0:
        # Prepare output filepath
        missing_out_fp = os.path.expanduser(out_dir.rstrip('/')) + "/" + "missing_snps_list_" + out_filename
        # Start with clean slate since we are appending
        if os.path.exists(missing_out_fp):
            os.remove(missing_out_fp)
        # Create and open file for writing
        cof_miss_vcf = open(missing_out_fp, 'a+')
        for i in vcf_header:
            cof_miss_vcf.write(i)
        for key in missing_snp_names.keys():
            cof_miss_vcf.write('\t'.join(missing_snp_names[key]) + '\n')


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
