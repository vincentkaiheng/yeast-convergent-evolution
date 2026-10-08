#!/bin/env python3
# -*- coding: UTF-8 -*-

import os
from glob import glob

def calculate_genome_stats(fna_file):
    total_length = 0
    gc_count = 0

    with open(fna_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq = line.strip().upper()
            total_length += len(seq)
            gc_count += seq.count('G') + seq.count('C')

    genome_size_mb = total_length / 1e6
    gc_content_percent = (gc_count / total_length * 100) if total_length > 0 else 0
    return genome_size_mb, gc_content_percent

def main():
    input_dir = './76_fna'
    output_file = 'genome_stat.tsv'
    
    fna_files = glob(os.path.join(input_dir, '*.fna'))
    
    with open(output_file, 'w') as out_f:
        out_f.write('Genome\tSize(Mb)\tGC(%)\n')
        for fna_file in sorted(fna_files):
            genome_name = os.path.basename(fna_file)
            size_mb, gc_percent = calculate_genome_stats(fna_file)
            out_f.write(f'{genome_name}\t{size_mb:.2f}\t{gc_percent:.2f}\n')

if __name__ == '__main__':
    main()
