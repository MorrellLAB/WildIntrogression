#!/usr/bin/env python3


"""
Written with Claude 3.7 Sonnnet through Copilot

Script to turn inferred introgression intervals into genetic distance estimates
"""

import argparse
import csv

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Calculate genetic distance for genomic intervals.')
    parser.add_argument('--csv', default='/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Figure_and_Table_Drafts/Figures_and_Tables_v4/Tables_and_Files_v3/Table S5 - wbdc_likely_introgressed_segments-breeding.csv',
                        help='Path to CSV file containing intervals')
    parser.add_argument('--bed', default='/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/recombination/merged_genetic_physical_map.bed',
                        help='Path to BED file containing genetic positions')
    parser.add_argument('--output', default='interval_genetic_distances.tsv',
                        help='Path to output file')
    return parser.parse_args()

def read_intervals(csv_file):
    """Read intervals from CSV file, skipping header and single-SNP intervals."""
    intervals = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        # Skip first two lines (header)
        next(reader)
        next(reader)

        for row in reader:
            # Skip single-SNP intervals (bp_length = 1)
            if int(row[8]) == 1:
                continue

            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            sample = row[3]
            mbp_length = float(row[9])

            intervals.append({
                'chromosome': chromosome,
                'start': start,
                'end': end,
                'sample': sample,
                'mbp_length': mbp_length
            })

    return intervals

def read_genetic_positions(bed_file):
    """Read genetic positions from BED file, organizing by chromosome."""
    positions_by_chr = {}

    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chromosome = fields[0]
            pos = int(fields[2])  # Use end position (1-based)
            genetic_pos = fields[4]

            # Skip positions without genetic map data
            if genetic_pos == '--':
                continue

            if chromosome not in positions_by_chr:
                positions_by_chr[chromosome] = []

            positions_by_chr[chromosome].append({
                'position': pos,
                'genetic_position': float(genetic_pos)
            })

    # Sort positions for each chromosome
    for chromosome in positions_by_chr:
        positions_by_chr[chromosome].sort(key=lambda x: x['position'])

    return positions_by_chr

def calculate_genetic_distances(intervals, positions_by_chr):
    """Calculate genetic distance for each interval."""
    results = []

    for interval in intervals:
        chromosome = interval['chromosome']
        start = interval['start']
        end = interval['end']
        
        # Skip if chromosome not in genetic positions
        if chromosome not in positions_by_chr:
            results.append({
                **interval,
                'genetic_distance': 'NA'
            })
            continue

        # Find markers within the interval
        markers_in_interval = []
        for pos in positions_by_chr[chromosome]:
            if start <= pos['position'] <= end:
                markers_in_interval.append(pos)

        # Calculate genetic distance (difference between min and max)
        if len(markers_in_interval) < 2:
            genetic_distance = 'NA'
        else:
            genetic_positions = [marker['genetic_position'] for marker in markers_in_interval]
            genetic_distance = max(genetic_positions) - min(genetic_positions)

        results.append({
            **interval,
            'genetic_distance': genetic_distance
        })

    return results

def write_output(results, output_file):
    """Write results to output file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("Chromosome\tStart\tEnd\tSample\tGenetic_Distance(cM)\tPhysical_Distance(Mbp)\n")
        
        # Write results
        for result in results:
            genetic_distance = result['genetic_distance']
            if genetic_distance != 'NA':
                genetic_distance = f"{genetic_distance:.2f}"
                
            f.write(f"{result['chromosome']}\t{result['start']}\t{result['end']}\t{result['sample']}\t{genetic_distance}\t{result['mbp_length']}\n")

def main():
    args = parse_arguments()
    
    # Read intervals from CSV
    print(f"Reading intervals from {args.csv}")
    intervals = read_intervals(args.csv)
    print(f"Read {len(intervals)} intervals (excluding single-SNP intervals)")
    
    # Read genetic positions from BED
    print(f"Reading genetic positions from {args.bed}")
    positions_by_chr = read_genetic_positions(args.bed)
    total_markers = sum(len(positions_by_chr[chr]) for chr in positions_by_chr)
    print(f"Read {total_markers} markers with genetic positions across {len(positions_by_chr)} chromosomes")
    
    # Calculate genetic distances
    print("Calculating genetic distances...")
    results = calculate_genetic_distances(intervals, positions_by_chr)

    # Write output
    write_output(results, args.output)
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()


