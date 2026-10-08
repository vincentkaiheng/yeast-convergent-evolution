#!/bin/env python3
# -*- coding: utf-8 -*-
"""Add per-species orthogroup counts to the genome summary."""

import csv
import io
import os
import sys

def read_tsv_skip_comments(path):
    with open(path, 'rt', encoding='utf-8') as f:
        lines = [l for l in f.readlines() if not l.lstrip().startswith('//')]
    text = ''.join(lines)
    return io.StringIO(text)

def build_mapping(orth_path, key_field='species', value_field='present_OG_count'):
    fh = read_tsv_skip_comments(orth_path)
    reader = csv.DictReader(fh, delimiter='\t')
    mapping = {}
    for r in reader:
        if r is None:
            continue
        key = r.get(key_field)
        if key is None:
            for k in r.keys():
                if k and k.lower() == key_field.lower():
                    key = r[k]
                    break
        val = r.get(value_field)
        if val is None:
            for k in r.keys():
                if k and k.lower() == value_field.lower():
                    val = r[k]
                    break
        if key is None:
            continue
        mapping[key] = val
    return mapping

def merge_counts(genome_path, mapping, genome_key='Genome', out_path=None, out_overwrite=False):
    fh = read_tsv_skip_comments(genome_path)
    reader = csv.DictReader(fh, delimiter='\t')
    rows = list(reader)
    if not rows:
        raise SystemExit("No rows found in genome table: {}".format(genome_path))

    header = reader.fieldnames[:] if reader.fieldnames else list(rows[0].keys())
    out_field = 'Gene family number'
    if out_field not in header:
        header.append(out_field)

    def lookup(name):
        if not name:
            return ''
        if name in mapping:
            return mapping[name]
        alt = name.replace(' ', '_')
        if alt in mapping:
            return mapping[alt]
        alt2 = name.replace('_', ' ')
        if alt2 in mapping:
            return mapping[alt2]
        for k in mapping:
            if k.lower() == name.lower():
                return mapping[k]
        return ''

    for r in rows:
        gname = r.get(genome_key) or r.get('species') or ''
        r[out_field] = lookup(gname)

    if out_path is None:
        base, ext = os.path.splitext(genome_path)
        out_path = base + '_with_gf_counts' + ext

    if os.path.exists(out_path) and not out_overwrite:
        raise SystemExit("Output already exists: {}".format(out_path))

    with open(out_path, 'wt', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=header, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return out_path

if __name__ == '__main__':
    # Example: ./genefamily_count/Orthogroups.presence_counts.tsv
    orth_path = '/home/yankh/76_converge/genomic_data/genefamily_count/Orthogroups.presence_counts.tsv'
    # Example: ./related_files/genome_stat_updated.tsv
    genome_path = '/home/yankh/76_converge/genomic_data/genome_stat_updated.tsv'

    if not os.path.exists(orth_path):
        raise SystemExit("Orthogroup file not found: {}".format(orth_path))
    if not os.path.exists(genome_path):
        raise SystemExit("Genome table not found: {}".format(genome_path))

    mapping = build_mapping(orth_path)
    out_path = merge_counts(genome_path, mapping)
    print("Output saved to:", out_path)
