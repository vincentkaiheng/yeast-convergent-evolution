#!/usr/bin/env python3

"""Build a term presence/absence matrix and exact intersection counts."""

import os
import glob
import argparse
import csv
from collections import defaultdict, OrderedDict

def read_go_summary(path):
    gos = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                continue
            gos.add(line)
    return gos

def main(input_dir, output_dir, write_intersections=True, include_go_lists=False):
    os.makedirs(output_dir, exist_ok=True)

    pattern = os.path.join(input_dir, "*_PFAM_summary.txt")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise SystemExit(f"No matching files found: {pattern}")

    taxa = []
    taxa_gos = OrderedDict()
    for p in paths:
        basename = os.path.basename(p)
        if not basename.endswith("_PFAM_summary.txt"):
            continue
        taxonomy = basename[:-len("_PFAM_summary.txt")]
        gos = read_go_summary(p)
        taxa.append(taxonomy)
        taxa_gos[taxonomy] = gos

    if not taxa:
        raise SystemExit("No PFAM summary files recognized.")

    all_gos = sorted(set().union(*taxa_gos.values()))

    out_table = os.path.join(output_dir, "PFAM_upset_table.tsv")
    with open(out_table, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        header = ["PFAM"] + [taxon + "-yeast" for taxon in taxa]
        writer.writerow(header)
        for go in all_gos:
            row = [go]
            for t in taxa:
                row.append("1" if go in taxa_gos[t] else "0")
            writer.writerow(row)

    print(f"Wrote UpSet table: {out_table}")

    if write_intersections:
        go_to_present = {}
        for go in all_gos:
            present = tuple(sorted([t for t in taxa if go in taxa_gos[t]]))
            go_to_present[go] = present

        combo_counts = defaultdict(set)
        for go, combo in go_to_present.items():
            combo_counts[combo].add(go)

        out_combo = os.path.join(output_dir, "PFAM_upset_intersections.tsv")
        with open(out_combo, 'w', newline='') as cf:
            writer = csv.writer(cf, delimiter='\t')
            writer.writerow(["combo", "size", "count"] + (["pfams"] if include_go_lists else []))

            sorted_combos = sorted(combo_counts.keys(), key=lambda c: (-len(c), ",".join(c)))
            for combo in sorted_combos:
                combo_label = "&".join(combo) if combo else "NONE"
                gos = sorted(combo_counts[combo])
                count = len(gos)
                if include_go_lists:
                    writer.writerow([combo_label, str(len(combo)), str(count), ",".join(gos)])
                else:
                    writer.writerow([combo_label, str(len(combo)), str(count)])

        print(f"Wrote exact intersection counts: {out_combo}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a term-by-lineage presence/absence matrix.")
    parser.add_argument("--input-dir", "-i",
                        default="./PFAMs",
                        help="Directory containing <lineage>_PFAM_summary.txt files")
    parser.add_argument("--output-dir", "-o",
                        default=".",
                        help="Output directory")
    parser.add_argument("--no-intersections", action="store_true",
                        help="Skip intersection counts")
    parser.add_argument("--include-pfam-lists", "--include-go-lists", dest="include_go_lists", action="store_true",
                        help="Include the terms in each intersection")
    args = parser.parse_args()

    main(args.input_dir, args.output_dir,
         write_intersections=not args.no_intersections,
         include_go_lists=args.include_go_lists)
