# Functional overlap across six yeast lineages

This directory compares gene-family changes across the six yeast lineages. Lost, expanded and contracted orthogroups are compared directly. Lost and novel families are also compared through their GO and PFAM annotations, allowing different OGs to contribute to the same functional category.

## Inputs from the preceding analyses

The third analysis directory supplies the PAPS lost and novel queries and the CAFE node-specific expansion and contraction lists. GO and PFAM terms come from the per-OG eggNOG-mapper files in `all_OGs_anno`, using the annotation workflow introduced in the first directory. Keep the same OrthoFinder OG identifiers throughout these steps.

The PAPS scripts use `HG_n -> OG(n - 1)`, with seven-digit OG numbers. This mapping assumes that the PAPS input retained the original OrthoFinder row order. The species tree and dated CAFE subtrees from the preceding analyses define the focal nodes; no new tree reconstruction is performed here.

Set the input paths near the start of each script and run it from its analysis subdirectory. Example tables and figures are included.

Use all six PAPS query files when interpreting intersections. The OG collectors warn and insert empty columns for missing groups; those columns do not establish biological absence.

## Lineage names

PAPS output filenames use the short lineage name, for example `Pucciniomycotina_OG_HG_GOs.tsv`. OG matrices use the same short names. GO/PFAM matrices retain the lost-analysis convention `<lineage>-yeast`; intersection combination labels use the short names.

The OG plots use the following node labels:

| Plot label | Yeast lineage | CAFE input node |
|---|---|---|
| Node 1 | Pucciniomycotina | node37 |
| Node 2 | Ustilaginomycotina | node45_basi |
| Node 3 | Agaricomycotina | node52 |
| Node 4 | Taphrinomycotina | node45_asco |
| Node 5 | Pezizomycotina | node82 |
| Node 6 | Saccharomycotina | node53 |

CAFE matrices retain their `<lineage>_yeast` column names. The `asco` and `basi` suffixes distinguish node 45 in the two separately fitted trees.

## 1. Shared lost orthogroups

Work in `PAPS_lost_intersections`.

1. Run `1_collect_upset_table.py`. It reads the PAPS lost-query gene lists, converts HG IDs to OG IDs, and writes the six-lineage presence/absence matrix in `OGs_upset_outputs`. It also writes separate tables for OGs occurring in exactly two through six lineage sets.
2. Run `plot_upset.py` to plot OGs shared by at least two lost-gene sets.
3. From `OGs_detail`, run `generate_OGs_detail.py` to list these shared OGs, their lineage combinations, and representative names, descriptions and PFAM strings. These detail columns summarize annotations; they are not the complete term sets used below.

The lost groups are absent from the focal yeast lineage and retained in its defined non-yeast sister comparison group, with the additional PAPS outgroup condition. Their selection is inherited from the third analysis directory.

## 2. GO and PFAM overlap among lost families

Continue in `PAPS_lost_intersections`.

1. Run `2_generate_go_table.py` and `2_generate_pfam_table.py`. For every selected OG, they collect all distinct GO terms from column 10 and PFAM annotations from column 21 of its eggNOG-mapper file. They write OG-level tables and one unique-term summary per lineage in `GOs` and `PFAMs`.
2. Run `3_collect_upset_table_go.py` and `3_collect_upset_table_pfam.py`. These convert the lineage summaries into binary matrices and exact intersection counts. Use `--include-go-lists` or `--include-pfam-lists`, respectively, to include the terms belonging to each combination.
3. Run `plot_GOs_upset.py` and `plot_PFAMs_upset.py` to plot terms shared by at least two of the six lineage sets.

A term is present in a lineage if at least one selected OG has that annotation. Repeated annotations within or across OGs do not increase its presence value. Missing annotations contribute no terms.

## 3. GO and PFAM overlap among novel families

Work in `PAPS_novel_intersections` and run the corresponding scripts in the same order: `1_collect_upset_table.py`, the two `2_generate_*` scripts, the two `3_collect_upset_table_*` scripts, and the two term-plotting scripts.

The input filters select the original **novel** category, not novel core: an OG must occur in at least half of the focal yeast species and be absent from every other sampled species. The thresholds are 3, 4, 3, 4, 6 and 2 for Pucciniomycotina, Ustilaginomycotina, Agaricomycotina, Taphrinomycotina, Saccharomycotina and Pezizomycotina, respectively.

These lineage-restricted OG sets are disjoint by definition. Their OG matrix records membership, and the exact-k OG tables are expected to be empty for k >= 2. The functional comparison instead asks whether different novel OGs carry the same GO or PFAM annotations.

## 4. Shared expansions and contractions

Work in `CAFE_intersections`. The supplied `node*.expanded` and `node*.contracted` lists represent the six nodes in the table above.

```sh
python 1_collect_OGs_upset_table.py --change expanded
python plot_upset.py --change expanded
python 1_collect_OGs_upset_table.py --change contracted
python plot_upset.py --change contracted
```

The collector writes separate OG matrices and nonempty exact-k intersection tables for each direction. It ignores the `FamilyID` input header. A present but empty list represents zero selected families; a missing node file is an input error. The plots retain all six node labels, including nodes with no shared families.

Run `generate_OGs_detail.py` from `OGs_expanded_detail` or `OGs_contracted_detail` to annotate OGs shared by at least two nodes. These scripts read the corresponding matrix and add representative names and descriptions. CAFE expansion and contraction are analyzed separately; these scripts do not impose an additional significance filter.

## 5. GO treemaps

Each PAPS `GOs` directory contains an R notebook, `treemap.ipynb`, with a REVIGO result for GO terms shared across all six lineages. Run it from that `GOs` directory with an R kernel and the `treemap` package to write `revigo_treemap.pdf`.

The REVIGO grouping step is external to these scripts. If the GO input changes, obtain an updated REVIGO result and replace the embedded data before rerunning the notebook. Rectangle area uses the supplied `uniqueness` value, and color identifies the representative term group.

## Reading the intersections

An exact intersection contains only items belonging to that particular combination of lineages. An exact-k table includes all combinations of k lineages and excludes items present in any additional lineage. The UpSet plots omit items found in only one lineage, so their set-size bars describe the shared subset rather than the full input lists.

## Software

Install the Python dependencies with `python -m pip install -r requirements.txt`. The version ranges avoid plotting incompatibilities between UpSetPlot 0.9 and newer pandas/NumPy releases. The GO treemap notebooks use R, an R Jupyter kernel and `treemap`.
