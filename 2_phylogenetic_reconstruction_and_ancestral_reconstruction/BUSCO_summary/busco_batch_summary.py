#!/bin/env python3

import os
from pathlib import Path
import pandas as pd

base_dir = Path(__file__).resolve().parent
busco_file = base_dir / "batch_summary.txt"
trait_file = (
    base_dir.parents[1] / "1_data_acquisition_and_functional_annotation"
    / "related_files" / "76_fungi_trait.txt"
)

output_dir = base_dir
output_file = os.path.join(output_dir, "busco_summary.txt")

os.makedirs(output_dir, exist_ok=True)

busco_df = pd.read_csv(busco_file, sep="\t")

trait_df = pd.read_csv(trait_file, sep="\t")

name_mapping = dict(zip(trait_df["file_names"], trait_df["names"]))

busco_df["Input_file"] = (
    busco_df["Input_file"]
    .astype(str)
    .str.replace(r"\.faa$", "", regex=True)
)

busco_df["Input_file"] = busco_df["Input_file"].map(name_mapping)

result_df = busco_df[
    ["Input_file", "Single", "Duplicated", "Fragmented"]
]

unmatched = result_df["Input_file"].isna().sum()
if unmatched > 0:
    print(f"Warning: {unmatched} entries could not be matched to species names.")

result_df.to_csv(output_file, sep="\t", index=False)

print(f"Output saved to: {output_file}")
