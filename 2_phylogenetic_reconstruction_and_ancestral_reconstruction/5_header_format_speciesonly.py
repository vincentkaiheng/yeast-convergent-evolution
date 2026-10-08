#!/bin/env python3
# -*- coding: UTF-8 -*-

import os

folder_path = './trimmed_sequences'
output_folder = './trimmed_sequences_sp_only'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for filename in os.listdir(folder_path):
    if filename.endswith('.faa'):
        file_path = os.path.join(folder_path, filename)
        
        with open(file_path, 'r') as file:
            lines = file.readlines()

        new_filename = os.path.join(output_folder, filename)
        with open(new_filename, 'w') as file:
            for line in lines:
                if line.startswith('>'):
                    # Headers must begin with the species identifier, followed by '|'.
                    line = line.split('|')[0] + '\n'
                file.write(line)

print("Updated FASTA files saved to:", output_folder)
