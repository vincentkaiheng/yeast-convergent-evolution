#!/bin/env python3
# -*- coding: UTF-8 -*-

from ete3 import Tree
import re
from pathlib import Path

tree_file = (Path(__file__).resolve().parents[3]
             / "2_phylogenetic_reconstruction_and_ancestral_reconstruction"
             / "MBASR" / "input" / "fungi_tree.nwk")
tree = Tree(str(tree_file))

tree_string = tree.write(format=0)

tree_string_cleaned = re.sub(r"(:\d+\.\d+|\.)", "", tree_string)
tree_string_cleaned = re.sub(r"(?<=\))\s*\d+\s*", "", tree_string_cleaned)

output_file = "./cleaned_fungi_tree.nwk"
with open(output_file, "w") as f:
    f.write(tree_string_cleaned)

print("Cleaned tree saved to:", output_file)
