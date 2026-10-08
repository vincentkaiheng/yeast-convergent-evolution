#!/usr/bin/env python3

import os

folder_path = '../before_RER/shared_orthologs'
output_folder = './shared_orthologs_sp_only'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for filename in os.listdir(folder_path):
    if filename.endswith('.fa'):
        file_path = os.path.join(folder_path, filename)
        
        with open(file_path, 'r') as file:
            lines = file.readlines()

        new_filename = os.path.join(output_folder, filename)
        with open(new_filename, 'w') as file:
            for line in lines:
                if line.startswith('>'):
                    line = line.split('|')[0] + '\n'
                file.write(line)

print("Species-only headers saved to:", output_folder)
