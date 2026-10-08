#!/bin/env python3
# -*- coding: UTF-8 -*-

import re

def clean_newick(newick_str):
    # Remove family counts and significance marks after each node ID.
    pattern = re.compile(r'(>)[^:]*')
    cleaned = pattern.sub(r'\1', newick_str)
    return cleaned

with open("tree2fig.nwk", "r", encoding="utf-8") as infile:
    tree_str = infile.read()

cleaned_tree = clean_newick(tree_str)

with open("tree2fig.nwk", "w", encoding="utf-8") as outfile:
    outfile.write(cleaned_tree)
