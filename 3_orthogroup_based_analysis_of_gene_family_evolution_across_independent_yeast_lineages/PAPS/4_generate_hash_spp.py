#!/bin/env python3
# -*- coding: UTF-8 -*-
"""Generate candidate PAPS taxonomy entries using NCBI Entrez."""
import os
import sys
import time
import argparse
from Bio import Entrez
from xml.etree import ElementTree as ET

SLOT_COUNT = 41  # Must match the nesting in PAPS.pl.
SLEEP_SEC = 1
RANK_ALIASES = {
    'superkingdom': ['superkingdom', 'superkingdom'],
    'kingdom': ['kingdom'],
    'phylum': ['phylum', 'division'],
    'subphylum': ['subphylum'],
    'class': ['class'],
    'subclass': ['subclass'],
    'order': ['order'],
    'suborder': ['suborder'],
    'infraorder': ['infraorder'],
    'family': ['family'],
    'genus': ['genus'],
}

SLOT_ASSIGNMENT = {
    0: 'superkingdom',
    1: 'node_a',
    2: 'node_b',
    3: 'node_c',
    4: 'phylum',
    5: 'phylum_b',
    6: 'phylum_c',
    7: 'phylum_d',
    8: 'phylum_e',
    9: 'phylum_f',
    10: 'phylum_g',
    11: 'phylum_h',
    12: 'phylum_i',
    13: 'phylum_j',
    14: 'phylum_k',
    15: 'subphylum_a',
    16: 'subphylum_b',
    17: 'subphylum_c',
    18: 'subphylum_d',
    19: 'superclass',
    20: 'class_a',
    21: 'class_b',
    22: 'class_c',
    23: 'class_d',
    24: 'class_e',
    25: 'class_f',
    26: 'class_g',
    27: 'subclass_a',
    28: 'subclass_b',
    29: 'subclass_c',
    30: 'order_a',
    31: 'order_b',
    32: 'order_c',
    33: 'order_d',
    34: 'suborder',
    35: 'infraorder',
    36: 'family_a',
    37: 'family_b',
    38: 'genus',
    39: 'species',
    40: "''"
}

def perlsinglequote(s):
    if s is None:
        return "''"
    s = s.replace("'", "\\'")
    return "'" + s + "'"

def read_ordered(ordered_path):
    with open(ordered_path, 'r', encoding='utf-8') as f:
        txt = f.read().strip()
    if not txt:
        return []
    parts = txt.split()
    return parts

def read_labelmap(mapping_path):
    m = {}
    name_from_abbr = {}
    with open(mapping_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                full, abbr = line.split("\t")[:2]
            else:
                parts = line.split()
                if len(parts) >= 2:
                    full = " ".join(parts[:-1])
                    abbr = parts[-1]
                else:
                    continue
            m[abbr] = full
            name_from_abbr[abbr] = full
    return name_from_abbr

def query_taxonomy(scientific_name, email):
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="taxonomy", term='"%s"[Scientific Name]' % scientific_name, retmode="xml")
        rec = Entrez.read(handle)
        handle.close()
        idlist = rec.get('IdList', [])
        if not idlist:
            time.sleep(0.2)
            handle = Entrez.esearch(db="taxonomy", term=scientific_name, retmode="xml")
            rec = Entrez.read(handle); handle.close()
            idlist = rec.get('IdList', [])
        if not idlist:
            return None
        taxid = idlist[0]
        time.sleep(SLEEP_SEC)
        fetch = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        xml = fetch.read()
        fetch.close()
        root = ET.fromstring(xml)
        rec0 = root.find('.//Taxon')
        if rec0 is None:
            return None
        found_name = rec0.findtext('ScientificName')
        lineage_ex = []
        for node in rec0.findall('./LineageEx/Taxon'):
            rn = node.findtext('Rank')
            sn = node.findtext('ScientificName')
            if rn is None: rn = ""
            if sn is None: sn = ""
            lineage_ex.append({'Rank': rn, 'ScientificName': sn})
        own_rank = rec0.findtext('Rank') or ''
        own_name = rec0.findtext('ScientificName') or ''
        return {'taxid': taxid, 'found_name': found_name, 'own_rank': own_rank, 'own_name': own_name, 'lineage_ex': lineage_ex}
    except Exception as e:
        return {'error': str(e)}

def build_slots_from_lineage(res, full_name, abbr):
    slots = ["''"] * SLOT_COUNT
    lineage_list = [ (t['Rank'].lower(), t['ScientificName']) for t in res.get('lineage_ex', []) ]
    lineage_names = [n for r,n in lineage_list if n]
    rk_map = { r: n for r,n in lineage_list }
    def get_rank(name_choices):
        for choice in name_choices:
            if choice in rk_map and rk_map[choice]:
                return rk_map[choice]
        return None

    sk = get_rank(['superkingdom', 'superkingdom', 'kingdom'])
    if sk:
        slots[0] = perlsinglequote(sk)
    else:
        if lineage_names:
            slots[0] = perlsinglequote(lineage_names[0])
        else:
            slots[0] = "''"

    nodes = []
    for r,nm in lineage_list:
        if nm and nm != slots[0].strip("'"):
            nodes.append(nm)
    nodes = [n for i,n in enumerate(nodes) if n not in nodes[:i]]
    for i in range(3):
        if i < len(nodes):
            slots[1 + i] = perlsinglequote(nodes[i])
        else:
            slots[1 + i] = "''"

    phyl = get_rank(['phylum','division'])
    if phyl:
        slots[4] = perlsinglequote(phyl)
    subp = get_rank(['subphylum'])
    if subp:
        slots[15] = perlsinglequote(subp)
    superc = get_rank(['superclass','subkingdom','kingdom'])
    if superc:
        slots[19] = perlsinglequote(superc)
    cls = get_rank(['class'])
    if cls:
        slots[20] = perlsinglequote(cls)
    ordn = get_rank(['order'])
    if ordn:
        slots[30] = perlsinglequote(ordn)
    fam = get_rank(['family'])
    if fam:
        slots[36] = perlsinglequote(fam)
    gen = get_rank(['genus'])
    if gen:
        slots[38] = perlsinglequote(gen)
    spname = res.get('own_name') or full_name
    if spname:
        sp_und = spname.replace(" ", "_")
        slots[39] = perlsinglequote(sp_und + "_(" + abbr + ")")
    else:
        slots[39] = perlsinglequote(abbr)

    return slots, lineage_names

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ordered", required=True, help="ordered_label.txt (one line, abbreviations separated by spaces)")
    parser.add_argument("--mapping", required=True, help="label.txt (tab-separated: full_name <TAB> abbr)")
    parser.add_argument("--out", default="hash_spp_generated.txt", help="Perl output file (default hash_spp_generated.txt)")
    parser.add_argument("--debug", default="hash_spp_debug.tsv", help="Debug TSV output")
    args = parser.parse_args()

    ordered = read_ordered(args.ordered)
    if not ordered:
        print("Error: ordered list is empty or not found.", file=sys.stderr)
        sys.exit(2)
    mapping = read_labelmap(args.mapping)

    email = os.environ.get('ENTREZ_EMAIL') or os.environ.get('ENTREZ_EMAIL'.lower())
    if not email:
        print("ERROR: set environment variable ENTREZ_EMAIL to a valid email (NCBI Entrez requirement).", file=sys.stderr)
        sys.exit(2)
    print("Using Entrez email:", email)
    out_lines = []
    debug_rows = []
    for idx, abbr in enumerate(ordered):
        full = mapping.get(abbr)
        if not full:
            print(f"[WARN] abbreviation {abbr} not found in mapping file; will generate empty placeholders.")
            slots = ["''"] * SLOT_COUNT
            slots[39] = perlsinglequote(abbr)
            slots[0] = perlsinglequote("Eukaryota")
            perl_line = "$spp " + "".join("{" + s + "}" for s in slots) + " = %d;" % idx
            out_lines.append(perl_line)
            debug_rows.append([idx,abbr, '', '', '', 'NOT_FOUND'] + slots)
            continue

        print(f"[{idx}] Querying NCBI taxonomy for: {full} ({abbr}) ...")
        res = query_taxonomy(full, email)
        if not res:
            print(f"[WARN] No taxonomy record found for {full} ({abbr}). Creating placeholders.")
            slots = ["''"] * SLOT_COUNT
            slots[39] = perlsinglequote(full.replace(" ", "_") + "_(" + abbr + ")")
            slots[0] = perlsinglequote("Eukaryota")
            perl_line = "$spp " + "".join("{" + s + "}" for s in slots) + " = %d;" % idx
            out_lines.append(perl_line)
            debug_rows.append([idx,abbr, full, '', 'NO_TAXON'] + slots)
            time.sleep(SLEEP_SEC)
            continue

        if isinstance(res, dict) and res.get('error'):
            print(f"[ERROR] Entrez error for {full}: {res['error']}")
            slots = ["''"] * SLOT_COUNT
            slots[39] = perlsinglequote(full.replace(" ", "_") + "_(" + abbr + ")")
            slots[0] = perlsinglequote("Eukaryota")
            perl_line = "$spp " + "".join("{" + s + "}" for s in slots) + " = %d;" % idx
            out_lines.append(perl_line)
            debug_rows.append([idx,abbr, full, '', 'ENTREZ_ERROR', res.get('error','')] + slots)
            time.sleep(SLEEP_SEC)
            continue

        slots, lineage_names = build_slots_from_lineage(res, full, abbr)
        perl_line = "$spp " + "".join("{" + s + "}" for s in slots) + " = %d;" % idx
        out_lines.append(perl_line)
        lineage_text = ";".join(lineage_names)
        debug_rows.append([idx,abbr, full, res.get('taxid',''), res.get('found_name',''), lineage_text] + slots)

        time.sleep(SLEEP_SEC)

    with open(args.out, 'w', encoding='utf-8') as fo:
        fo.write("# Generated by 4_generate_hash_spp.py\n")
        fo.write("# Review taxonomy and restore custom lineage/sister groups before updating hash_spp.\n\n")
        for ln in out_lines:
            fo.write(ln + "\n")
    print("Perl lines written to", args.out)

    with open(args.debug, 'w', encoding='utf-8') as fd:
        headers = ["index","abbr","full_name","taxid","found_name","lineage_text"] + [f"slot{i}" for i in range(SLOT_COUNT)]
        fd.write("\t".join(headers) + "\n")
        for row in debug_rows:
            fd.write("\t".join(str(x) for x in row) + "\n")
    print("Debug TSV written to", args.debug)
    print("Done. Inspect debug TSV, manually review any entries with NOT_FOUND / NO_TAXON / ENTREZ_ERROR. Preserve custom lineage and sister-group definitions when updating PAPS.pl.")

if __name__ == "__main__":
    main()
