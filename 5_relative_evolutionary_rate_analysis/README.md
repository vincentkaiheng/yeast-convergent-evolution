# Relative evolutionary rates in yeast and non-yeast fungi

This workflow filters OrthoFinder orthogroups, prunes paralogs, and tests whether protein evolutionary rates are associated with the yeast growth form. Positive associations identify relatively accelerated ortholog sets; negative associations identify relatively decelerated sets. Example inputs and selected results are included.

## Inputs from earlier analyses

The starting gene-count matrix and protein sequences come from the OrthoFinder analysis in directory 3. Protein headers use `species_id|sequence_id`.

The shared species table is `../1_data_acquisition_and_functional_annotation/related_files/76_fungi_trait.txt`. Its `numbering` and `names` columns connect numeric sequence identifiers to species names. The 43 identifiers in `before_RER/foreground_list.txt` match the species with `trait = 1` in that table.

`RER/related_files/fungi_tree.tre` is the species-tree input, prepared from the BUSCO tree reconstructed in directory 2. It uses numeric tip labels, matching the sequence headers before renaming. Functional annotations come from the eggNOG-mapper results generated earlier in the workflow.

## 1. Filter and prune orthogroups

Run the scripts in order from `before_RER/`.

1. **`1_extract_12_OGs.py`** retains OGs present in more than 12 species, meaning at least 13, and copies their protein sequences to `OGs_13filtered`. The OrthoFinder `Total` column is excluded from the species count.
2. **`2_pre_pp_pipeline.sh`** aligns each OG with MAFFT, trims it with trimAl, and estimates a gene tree with FastTreeMP. The settings are MAFFT `--localpair --maxiterate 1000`, trimAl `-automated1`, and FastTreeMP `-lg -gamma`. Alignments and matching trees in `trimmed_OGs` are passed to the next step.
3. **`3_phylopypruner.sh`** uses PhyloPyPruner to remove paralogs with the maximum-inclusion method (`--prune MI`), a minimum sequence length of 50, minimum support of 0.50, and at least 13 taxa.
4. **`4_OG_seperator.py`** separates the pruned alignments into `yeast_only_orthologs`, `nonyeast_only_orthologs`, and `shared_orthologs`. Only the last group, containing both yeast and non-yeast species, proceeds to RER analysis.

One original OG may yield several pruned ortholog sets. Suffixes such as `_pruned_1` and `_pruned_2` identify separate sets and must be retained throughout the analysis.

## 2. Estimate relative evolutionary rates

Run the following steps from `RER/`.

1. **`1_header_format_sponly.py`** reads `before_RER/shared_orthologs` and keeps only the species identifier in each sequence header. The resulting alignments are written to `shared_orthologs_sp_only`.
2. **`2_phangorn.R`** estimates branch lengths for each alignment using the supplied species-tree topology and writes `related_files/gene_trees.tre`. Run it with `Rscript 2_phangorn.R`.
3. **`3_update_gene_tree.py`** replaces numeric tree-tip identifiers with species names using the shared trait table. Ortholog-set identifiers are preserved in `related_files/new_gene_trees.tre`.
4. **`4_count_minTreesAll.py`** counts alignments with 76 sequences. After single-copy pruning, these are the alignments containing all study species. The supplied gene-tree file has 428 such trees; this value is used for `minTreesAll` in the next step. Recheck it if the input alignments change.
5. **`5_RER.ipynb`** reads the gene trees, estimates relative rates with square-root transformation, weighting and scaling, and correlates them with the binary foreground using `min.sp = 10`, `min.pos = 2`, and `weighted = "auto"`.

The notebook defines the 43 yeast foreground species and their internal foreground branches. Its seven clade groups include separate Exophiala and Aureobasidium groups within Pezizomycotina; Acaromyces is a separate foreground tip. Glomus and Umbelopsis are used to root the foreground-tree display.

The output `cor_all.txt` contains 7,281 tested ortholog sets with usable statistics. Candidate sets are selected using **Rho > 0.35 and raw P < 0.05** for acceleration, or **Rho < -0.35 and raw P < 0.05** for deceleration. These give 277 accelerated and 161 decelerated sets in `cor_pos.txt` and `cor_neg.txt`.

`Rho` is the association coefficient returned by the binary-phenotype analysis, `N` is the number of branches included, and `stat` is `sign(Rho) × -log10(P)`.

## 3. Collect annotations and plot the results

1. **`6_extract_gene_files2anno.py`** uses the selected ortholog-set IDs to recover their original, unaligned protein sequences. It matches the full sequence IDs in the pruned alignments to the original OG files. Set `TABLE_PATH` and `OUT_DIR` for the positive or negative candidate set.
2. **`7_extract_annotation_with_dic.py`** retrieves eggNOG-mapper annotations for those sequences from the earlier annotation tables. Set the input and output directories for `pos2anno`/`pos_anno` or `neg2anno`/`neg_anno`.
3. **`8_collect_anno_with_pfam.py`** summarizes Description, Preferred_name and the reported PFAM combination for each tested set. It ignores annotation headers and missing values, chooses the most frequent informative value in each field, and uses the first encountered value to resolve ties. A dash marks an unavailable annotation. Run it for both directions by changing the three paths in `main()`.
4. **`RER_plot.ipynb`** plots Rho against the raw P value. The widgets adjust highlighted labels and save the plot and label positions. Its thresholds match the candidate lists above.

## Software

The workflow uses Python with pandas and Biopython; MAFFT, trimAl, FastTreeMP and PhyloPyPruner; R with RERconverge and its dependencies; and Jupyter with an R kernel. The plotting notebook additionally uses matplotlib and ipywidgets. The shell alignment pipeline requires Bash 4 or later. Adjust input paths and parallel job counts for the analysis environment.

The [RERconverge walkthrough](https://github.com/nclark-lab/RERconverge/blob/master/vignettes/BioinformaticsRERconvergeSupp.Rmd) describes the rate calculation and association statistics.
