# Phylogenetic reconstruction and ancestral state reconstruction

This workflow extracts single-copy BUSCO proteins, filters and aligns loci, builds a species tree, summarizes BUSCO completeness, and reconstructs ancestral growth states with MBASR. Included files illustrate the analysis; the full BUSCO and alignment datasets are not included.

## Connection to the previous directory

The shared species table is `../1_data_acquisition_and_functional_annotation/related_files/76_fungi_trait.txt`. It provides species names, numeric identifiers, protein-file names, and growth states (`1` = general yeast, `0` = multicellular fungus).

`BUSCO_summary/busco_batch_summary.py` uses this table to map protein-file names to species names. The single-trait MBASR example reads the same table and writes its `names` and `trait` columns to `MBASR/input/fungi_trait.txt`, without a header, as required by MBASR. Keep the two analysis directories alongside each other.

## Species-tree reconstruction

Run the numbered scripts from this directory after obtaining per-species BUSCO results with the `fungi_odb12.2` lineage dataset.

1. **`1_extract_single_copy_busco.py`** collects single-copy protein sequences from each species' BUSCO output and combines sequences with the same BUSCO identifier into one FASTA file per locus. Original sequence headers are retained.
2. **`2_single_copy_filter.py`** retains loci with at least 61 sequences, approximately 80% of the 76 species. The threshold can be changed with `--min-sequences`. This counts FASTA records and assumes one sequence per species in each locus.
3. **`3_parallel_mafft.sh`** aligns each retained locus with MAFFT using `--localpair --maxiterate 1000`.
4. **`4_parallel_trimal.sh`** trims the alignments with trimAl using `-automated1`.
5. **`5_header_format_speciesonly.py`** keeps the part of each FASTA header before the first `|`, leaving the species identifier for tree construction. Input headers must already follow this convention.
6. **`6_iqtree.sh`** passes the alignment directory to IQ-TREE for concatenated tree inference, with ModelFinder (`-m MFP`), 30 threads, and 1,000 ultrafast bootstrap replicates (`-B 1000`). The tree and run outputs are in `final_treefile`.

The IQ-TREE tips use numeric IDs. Map them to species names with `76_fungi_trait.txt` and root the tree before ancestral reconstruction. The prepared tree is `MBASR/input/fungi_tree.nwk`.

## BUSCO summary

`BUSCO_summary/busco_batch_summary.py` reads `batch_summary.txt`, removes the `.faa` suffix from input-file names, maps them to species names, and retains the `Single`, `Duplicated`, and `Fragmented` percentages. It writes `busco_summary.txt` in the same folder. This summary can be produced independently of the alignment and tree steps.

## Ancestral state reconstruction

Run `MBASR/1.example.commands.single.trait.txt` in R from this directory. It prepares the trait file and runs MBASR with `character.type = "ordered"` and `n.samples = 5000`.

MBASR results and plots are saved in `MBASR/output`. Edit `MBASR/input/plot.settings.txt` and use `replot()` to adjust the figure.

## Software

The Python scripts use Biopython and pandas. Tree construction requires BUSCO results, MAFFT, trimAl, and IQ-TREE; adjust parallel job counts to the available resources. MBASR uses R with ape, phytools, and pdftools, and checks for MrBayes 3.2.7 or 3.2.7a. Its loader installs missing R packages. Place the MrBayes executable in `MBASR/other/mb` as `mb` (macOS/Linux) or `mb.exe` (Windows).

## Toolkit sources

MBASR is from [Steven Heritage’s releases](https://github.com/stevenheritage/MBASR/releases), associated with [this preprint](https://www.biorxiv.org/content/10.1101/2021.01.10.426107v1.full). Its upstream CC0 license is included in `MBASR/LICENSE`. Obtain MrBayes separately from the [official download page](https://nbisweden.github.io/MrBayes/download.html); see `MBASR/other/mb/README.md`.
