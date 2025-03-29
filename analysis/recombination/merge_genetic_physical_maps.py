#!/usr/bin/env python3
"""
Written with Claude 3.7 Sonnnet through Copilot

Script to merge genetic and physical map data from a genetic positions file
and a VCF file, based on common marker IDs.

Can output in the following formats:
1. Tab-separated file containing:
   - Marker ID
   - Chromosome
   - Genetic position
   - Physical position
   - Reference allele
   - Alternate allele

2. BED format file containing:
   - Chromosome (with 'chr' prefix)
   - Start position (physical position - 1)
   - End position (physical position)
   - Marker ID
   - Genetic position (or '--' if not available)

Usage:
  python3 merge_genetic_physical_maps.py --format [bed|tsv|both]
"""

import sys
import argparse

# Input file paths
GENETIC_POS_FILE = "bopa_and_9k_genetic_pos.txt"
VCF_FILE = "bopa_idt90_noRescuedSNPs.vcf"

# Output file paths
OUTPUT_TSV = "merged_genetic_physical_map.tsv"
OUTPUT_BED = "merged_genetic_physical_map.bed"


def load_genetic_positions(file_path):
    """Load genetic positions from the specified file."""
    genetic_data = {}
    chr_counts = {}

    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Skip empty lines
                if line.strip() == "":
                    continue

                # Parse the line
                parts = line.strip().split()
                if len(parts) < 3:
                    print(f"Warning: Skipping malformed line in genetic positions file: {line.strip()}")
                    continue

                marker_id = parts[0]
                chromosome = parts[1]
                genetic_pos = float(parts[2])

                # Store the data
                genetic_data[marker_id] = {
                    'chromosome': chromosome,
                    'genetic_pos': genetic_pos
                }

                # Count chromosomes
                chr_counts[chromosome] = chr_counts.get(chromosome, 0) + 1

        print(f"Loaded {len(genetic_data)} markers from genetic positions file")
        print("Chromosomes in genetic data:", ", ".join(f"{chr}({count})" for chr, count in sorted(chr_counts.items())))

        return genetic_data

    except Exception as e:
        print(f"Error loading genetic positions file: {e}")
        sys.exit(1)

def load_vcf_data(file_path):
    """Load physical positions and alleles from VCF file."""
    vcf_data = {}
    chr_counts = {}

    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue

                # Parse the line
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    print(f"Warning: Skipping malformed line in VCF file: {line.strip()}")
                    continue

                chromosome = parts[0]
                physical_pos = int(parts[1])
                marker_id = parts[2]
                ref_allele = parts[3]
                alt_allele = parts[4]

                # Remove 'chr' prefix if present for consistency
                clean_chromosome = chromosome.replace('chr', '')

                # Store the data
                vcf_data[marker_id] = {
                    'chromosome': clean_chromosome,
                    'physical_pos': physical_pos,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                    'original_chromosome': chromosome
                }

                # Count chromosomes
                chr_counts[clean_chromosome] = chr_counts.get(clean_chromosome, 0) + 1

        print(f"Loaded {len(vcf_data)} markers from VCF file")
        print("Chromosomes in VCF data:", ", ".join(f"{chr}({count})" for chr, count in sorted(chr_counts.items())))

        return vcf_data

    except Exception as e:
        print(f"Error loading VCF file: {e}")
        sys.exit(1)

def merge_data(genetic_data, vcf_data):
    """Merge genetic and physical data based on ALL markers from VCF file."""
    merged_data = []
    chr_match_count = 0
    chr_mismatch_count = 0
    markers_with_genetic_pos = 0
    markers_without_genetic_pos = 0

    # Find common markers
    common_markers = set(genetic_data.keys()) & set(vcf_data.keys())
    print(f"Found {len(common_markers)} common markers between the two datasets")
    print(f"Total markers in VCF file: {len(vcf_data)}")

    # Process all VCF markers
    for marker_id, vcf_info in vcf_data.items():
        # Check if marker has genetic position data
        has_genetic_data = marker_id in genetic_data

        if has_genetic_data:
            genetic_info = genetic_data[marker_id]

            # Check chromosome consistency
            chr_match = genetic_info['chromosome'] == vcf_info['chromosome']
            if chr_match:
                chr_match_count += 1
            else:
                chr_mismatch_count += 1
                print(f"Chromosome mismatch for marker {marker_id}: "
                      f"Genetic={genetic_info['chromosome']}, "
                      f"VCF={vcf_info['chromosome']} (original={vcf_info['original_chromosome']})")

            # Create merged record with genetic position
            merged_record = {
                'marker_id': marker_id,
                'chromosome': genetic_info['chromosome'] if chr_match else vcf_info['chromosome'],
                'genetic_pos': genetic_info['genetic_pos'],
                'has_genetic_pos': True,
                'physical_pos': vcf_info['physical_pos'],
                'ref_allele': vcf_info['ref_allele'],
                'alt_allele': vcf_info['alt_allele'],
                'chr_match': chr_match
            }
            markers_with_genetic_pos += 1
        else:
            # Create record without genetic position
            merged_record = {
                'marker_id': marker_id,
                'chromosome': vcf_info['chromosome'],
                'genetic_pos': None,
                'has_genetic_pos': False,
                'physical_pos': vcf_info['physical_pos'],
                'ref_allele': vcf_info['ref_allele'],
                'alt_allele': vcf_info['alt_allele'],
                'chr_match': True  # No comparison possible
            }
            markers_without_genetic_pos += 1

        merged_data.append(merged_record)

    print(f"Chromosome matches: {chr_match_count}")
    print(f"Chromosome mismatches: {chr_mismatch_count}")
    print(f"Markers with genetic positions: {markers_with_genetic_pos}")
    print(f"Markers without genetic positions: {markers_without_genetic_pos}")
    
    return merged_data

def write_tsv_output(merged_data, output_file):
    """Write merged data to tab-delimited output file."""
    try:
        with open(output_file, 'w') as f:
            # Write header
            header = "Marker_ID\tChromosome\tGenetic_Position\tPhysical_Position\tReference_Allele\tAlternate_Allele\n"
            f.write(header)

            # Sort by chromosome and then by genetic position
            # Filter out records without genetic positions for TSV output
            filtered_data = [record for record in merged_data if record['has_genetic_pos']]

            # Sort by chromosome and then by genetic position
            sorted_data = sorted(filtered_data, 
                                key=lambda x: (x['chromosome'], x['genetic_pos']))

            # Write data
            for record in sorted_data:
                line = (f"{record['marker_id']}\t{record['chromosome']}\t"
                        f"{record['genetic_pos']}\t{record['physical_pos']}\t"
                        f"{record['ref_allele']}\t{record['alt_allele']}\n")
                f.write(line)

        print(f"Wrote {len(sorted_data)} merged markers with genetic positions to {output_file}")

    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

def write_bed_output(merged_data, output_file):
    """Write merged data to BED format output file."""
    try:
        with open(output_file, 'w') as f:
            # Sort by chromosome and then by position
            sorted_data = sorted(merged_data, 
                                key=lambda x: (x['chromosome'], x['physical_pos']))

            # Write data
            for record in sorted_data:
                # BED format: chrom, start, end, name, score, strand, etc.
                chrom = f"chr{record['chromosome']}"
                start = record['physical_pos'] - 1  # BED is 0-based, start is inclusive
                end = record['physical_pos']  # BED end is exclusive
                name = record['marker_id']
 
                # Use raw genetic position or '--' if not available
                if record['has_genetic_pos']:
                    score = str(record['genetic_pos'])
                else:
                    score = '--'

                line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\n"
                f.write(line)

        print(f"Wrote {len(merged_data)} merged markers to {output_file}")

    except Exception as e:
        print(f"Error writing BED output file: {e}")
        sys.exit(1)

def main():
    """Main function to execute the merge workflow."""
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description='Merge genetic and physical map data from genetic positions file and VCF file.')
    parser.add_argument('--format', choices=['bed', 'tsv', 'both'], default='bed',
                        help='Output format: tab-delimited (tsv), BED format (bed), or both')
    parser.add_argument('--genetic', default=GENETIC_POS_FILE,
                        help=f'Path to genetic positions file (default: {GENETIC_POS_FILE})')
    parser.add_argument('--vcf', default=VCF_FILE,
                        help=f'Path to VCF file (default: {VCF_FILE})')
    parser.add_argument('--out-tsv', default=OUTPUT_TSV,
                        help=f'Path to tab-delimited output file (default: {OUTPUT_TSV})')
    parser.add_argument('--out-bed', default=OUTPUT_BED,
                        help=f'Path to BED output file (default: {OUTPUT_BED})')

    args = parser.parse_args()

    print(f"Loading genetic positions from {args.genetic}...")
    genetic_data = load_genetic_positions(args.genetic)

    print(f"\nLoading physical positions and alleles from {args.vcf}...")
    vcf_data = load_vcf_data(args.vcf)

    print("\nMerging data based on common marker IDs...")
    merged_data = merge_data(genetic_data, vcf_data)
 
    # Output the data in the specified format(s)
    if args.format in ['tsv', 'both']:
        print(f"\nWriting merged data with genetic positions to tab-delimited file {args.out_tsv}...")
        write_tsv_output(merged_data, args.out_tsv)

    if args.format in ['bed', 'both']:
        print(f"\nWriting ALL markers to BED format file {args.out_bed}...")
        write_bed_output(merged_data, args.out_bed)

    print("\nDone!")

if __name__ == "__main__":
    main()

