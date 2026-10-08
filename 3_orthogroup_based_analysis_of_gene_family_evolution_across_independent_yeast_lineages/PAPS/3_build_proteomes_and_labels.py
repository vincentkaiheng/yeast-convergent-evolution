#!/bin/env python3
# -*- coding: UTF-8 -*-

import os

order_file = "./order.txt"
label_file = "./label.txt"
faa_dir = "./faa"
output_faa = "./run_PAPS(modified_from_Jordi_Paps)/Input/all_proteomes.faa"
output_label = "./ordered_label.txt"

with open(order_file, "r") as f:
    order_list = [line.strip() for line in f if line.strip()]

label_dict = {}
with open(label_file, "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            species, abbr = parts[0], parts[1]
            label_dict[species] = abbr

os.makedirs(os.path.dirname(output_faa), exist_ok=True)

with open(output_faa, "w") as out_f:
    for species in order_list:
        faa_path = os.path.join(faa_dir, f"{species}.faa")
        if not os.path.exists(faa_path):
            raise FileNotFoundError(f"Missing FAA file: {faa_path}")
        with open(faa_path, "r") as sf:
            out_f.write(sf.read())

abbr_list = []
for species in order_list:
    if species not in label_dict:
        raise KeyError(f"Species not found in label file: {species}")
    abbr_list.append(label_dict[species])

with open(output_label, "w") as f:
    f.write(" ".join(abbr_list) + "\n")
