#!/usr/bin/env bash
# Retain pruned ortholog sets with at least 13 taxa.
directory=./trimmed_OGs
phylopypruner --dir "${directory}" --min-len 50 \
--min-support 0.50 --prune MI --min-taxa 13
