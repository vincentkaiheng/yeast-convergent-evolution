#!/bin/env python3
# -*- coding: UTF-8 -*-


import re

def scale_newick_tree(input_file, output_file, scale_factor):
    with open(input_file, 'r') as file:
        tree_str = file.read().strip()
    
    def scale_match(match):
        return f":{float(match.group(1)) * scale_factor}"
    
    scaled_tree = re.sub(r":([0-9]+\.?[0-9]*)", scale_match, tree_str)
    
    with open(output_file, 'w') as file:
        file.write(scaled_tree)

# Scale each prepared subtree; change asco to basi for the second run.
input_file = "./asco/tree_asco.txt"
output_file = "./asco/tree_asco_scaled.txt"

scale_newick_tree(input_file, output_file, 100)
