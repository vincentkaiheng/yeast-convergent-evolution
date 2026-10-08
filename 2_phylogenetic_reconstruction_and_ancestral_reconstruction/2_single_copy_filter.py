#!/bin/env python3
# -*- coding: UTF-8 -*-

import os
import sys

def count_sequences_in_file(file_path):
    """Count FASTA records."""
    try:
        with open(file_path, 'r') as f:
            count = 0
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return 0

def filter_fasta_files(input_dir, output_dir, min_sequences=61, quiet=False):
    """Retain loci with at least min_sequences FASTA records."""
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.faa')]
    total_files = len(fasta_files)

    if not quiet:
        print(f"Input directory: {input_dir}")
        print(f"Total files: {total_files}")
        print(f"Minimum sequences: {min_sequences}")
        print("-" * 50)

    filtered_count = 0
    discarded_count = 0
    sequence_counts = []

    for filename in fasta_files:
        file_path = os.path.join(input_dir, filename)
        seq_count = count_sequences_in_file(file_path)
        sequence_counts.append(seq_count)

        if seq_count >= min_sequences:
            output_path = os.path.join(output_dir, filename)

            with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
                outfile.write(infile.read())

            filtered_count += 1
            if not quiet:
                print(f"Keep {filename}: {seq_count} sequences")
        else:
            discarded_count += 1
            if not quiet:
                print(f"Skip {filename}: {seq_count} sequences (<{min_sequences})")

    print("\n" + "=" * 50)
    print("Filtering complete.")
    print("=" * 50)
    print(f"Input files: {total_files}")
    print(f"Retained files: {filtered_count} (>={min_sequences} sequences)")
    print(f"Discarded files: {discarded_count} (<{min_sequences} sequences)")

    if sequence_counts:
        print(f"\nSequence counts:")
        print(f"  Maximum: {max(sequence_counts)} sequences")
        print(f"  Minimum: {min(sequence_counts)} sequences")
        print(f"  Mean: {sum(sequence_counts)/len(sequence_counts):.1f} sequences")

        print(f"\nSequence-count distribution:")
        ranges = [(0, 20), (21, 40), (41, 60), (61, 80), (81, 100), (101, 120), (121, float('inf'))]
        for r_start, r_end in ranges:
            if r_end == float('inf'):
                count = sum(1 for c in sequence_counts if c >= r_start)
                print(f"  ≥{r_start}: {count} files")
            else:
                count = sum(1 for c in sequence_counts if r_start <= c <= r_end)
                print(f"  {r_start}-{r_end}: {count} files")

    print(f"\nRetained files saved to: {output_dir}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Retain FASTA files with at least the specified number of sequences.')
    parser.add_argument('-i', '--input', default='./combined_single_copy_sequences',
                       help='Input directory (default: ./combined_single_copy_sequences)')
    parser.add_argument('-o', '--output', default='./filtered_single_copy_sequences',
                       help='Output directory (default: ./filtered_single_copy_sequences)')
    parser.add_argument('-m', '--min-sequences', type=int, default=61,
                       help='Minimum sequences (default: 61, approximately 80% of 76 species)')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Print the summary only')

    args = parser.parse_args()

    filter_fasta_files(args.input, args.output, args.min_sequences, args.quiet)
