# Gene-family evolution across independent yeast lineages

OrthoFinder supplies the orthogroups for two analyses: PAPS screens their presence and absence across yeast lineages, while CAFE estimates family expansion and contraction on dated trees. CAFE follows time-tree reconstruction; PAPS can run independently of dating.

## Connections to the preceding analyses

- Protein FASTA files come from the dataset used in `1_data_acquisition_and_functional_annotation`. The shared `76_fungi_trait.txt` in that directory's `related_files` supplies species names, numeric IDs, file names, and yeast/non-yeast assignments.
- Dating uses a concatenation of the BUSCO alignments in `trimmed_sequences_sp_only` from `2_phylogenetic_reconstruction_and_ancestral_reconstruction`. Its species tree is prepared by removing unnecessary annotations and marking calibration nodes.
- The OrthoFinder gene-count table also supplies the gene-family counts summarized in the first directory.

The included files illustrate the workflow. Some large PAPS inputs, including the combined proteome, are omitted. Set the protein and orthogroup-sequence paths for the full dataset before running.

## OrthoFinder

`orthofinder.sh` runs OrthoFinder on the protein files with DIAMOND and the MSA workflow. `orthofinder_result` holds the orthogroup membership and gene-count tables used below. Run the script where `76_faa` contains the full protein dataset, or adjust its input path.

## PAPS

Run the preparation scripts from `PAPS`:

1. **`1_extract_rename_fasta.py`** replaces numeric species prefixes in protein headers with short labels, writes renamed protein files to `faa`, and creates `label.txt`. The label-file keys are filename stems without `.faa`; these must match `order.txt`.
2. **`2_orthofinder2mcl-like.py`** converts the OrthoFinder membership table into one space-separated list of protein IDs per orthogroup, using the same species labels. The converted file is also named `Orthogroups.tsv`, but is no longer the original tabular OrthoFinder format.
3. **`3_build_proteomes_and_labels.py`** combines the renamed proteins in `order.txt` order into `run_PAPS(modified_from_Jordi_Paps)/Input/all_proteomes.faa` and writes `ordered_label.txt`.
4. **`4_generate_hash_spp.py`** optionally retrieves NCBI taxonomy to draft the species entries for `PAPS.pl`. Run it with `--ordered ordered_label.txt --mapping label.txt` after setting `ENTREZ_EMAIL`. Review the generated taxonomy before use. Species groups are defined in `PAPS.pl`.
5. Place the converted `Orthogroups.tsv` in the PAPS runner's `Input` directory. From that directory, run `perl MCL_row_counter.pl Orthogroups.tsv` to create `gene_numbers_parsed.out`. Its species-label order must match `ordered_label.txt` and the column indices in `PAPS.pl`.
6. From `run_PAPS(modified_from_Jordi_Paps)`, run `perl Create_DBs.pl`, then `perl PAPS.pl` and enter the presence/absence queries below. Results are written to `Output`.

The PAPS Perl scripts are adapted from [JLWei7/animal_terrestrialisation](https://github.com/JLWei7/animal_terrestrialisation). PAPS assigns row-based IDs such as `HG_1`, so preserve the converted table's row order when matching results back to OrthoFinder orthogroups.

### Novel and novel core

Both categories refer specifically to the **yeast species within the focal subphylum**, compared with **all other sampled species**, including yeasts from other subphyla.

- **Novel:** present in at least half of the focal yeast species, rounding up, and absent from every species outside that focal yeast group.
- **Novel core:** present in at least `N - 1` of the `N` focal yeast species, and absent from every species outside that group. It permits absence from one focal species; it does not require presence in all of them.

Use `<lineage>_yeast-atleast<K> non_<lineage>_yeast-absent`, replacing the lineage and threshold with the entries below. Counts are orthogroups (PAPS homology groups), not individual genes.

| Focal yeast lineage | N | Novel K | Novel groups | Novel core K | Novel core groups |
|---|---:|---:|---:|---:|---:|
| Pucciniomycotina | 5 | 3 | 785 | 4 | 565 |
| Ustilaginomycotina | 8 | 4 | 336 | 7 | 17 |
| Agaricomycotina | 6 | 3 | 1,062 | 5 | 620 |
| Taphrinomycotina | 8 | 4 | 463 | 7 | 2 |
| Saccharomycotina | 12 | 6 | 226 | 11 | 26 |
| Pezizomycotina | 4 | 2 | 683 | 3 | 316 |

For example, `Pucciniomycotina_yeast-atleast3 non_Pucciniomycotina_yeast-absent` selects the 785 novel groups. These are lineage-restricted presence/absence candidates; the query itself does not establish how a family originated.

### Lost

Lost groups are absent from every yeast species in the focal lineage but retained in its defined sister comparison group. **Yeast species within the sister lineage are excluded from that comparison group.** The supplied queries require presence in every member of the corresponding `sister` group (`present` in PAPS), plus presence in at least one remaining sampled species (`outg-atleast1`). Here `outg` is the set left after excluding the focal and sister groups.

| Focal yeast lineage | Query | Lost groups |
|---|---|---:|
| Pucciniomycotina | `Pucciniomycotina_yeast-absent sister1-present outg-atleast1` | 228 |
| Ustilaginomycotina | `Ustilaginomycotina_yeast-absent sister2-present outg-atleast1` | 14 |
| Agaricomycotina | `Agaricomycotina_yeast-absent sister3-present outg-atleast1` | 328 |
| Taphrinomycotina | `Taphrinomycotina_yeast-absent sister4-present outg-atleast1` | 322 |
| Saccharomycotina | `Saccharomycotina_yeast-absent sister5-present outg-atleast1` | 210 |
| Pezizomycotina | `Pezizomycotina_yeast-absent sister6-present outg-atleast1` | 428 |

Each query produces corresponding files for species counts, gene IDs, sequence-header annotations, and taxon names. The query and number of selected groups appear in the filenames. Novel/core comparisons use all species outside the focal yeast group; lost comparisons use the specified non-yeast sister group and the additional outgroup condition.

## Time-tree reconstruction

Work in `CAFE/timetree`. `input.fa` contains the concatenated BUSCO protein alignment. `input.tre` contains the prepared species topology and calibration labels. Sequence concatenation and calibration placement precede these runs and are not automated by the supplied scripts. `clean_tree.py` assists with removing branch-length/support information from the earlier species tree; calibration labels still need to be added to the dating input.

1. Set `usedata = 3` in `mcmctree.ctl` and run `mcmctree mcmctree.ctl`. This generates `out.BV` and intermediate files for the likelihood approximation.
2. Use the supplied `wag.dat` matrix and edit the generated `tmp0001.ctl` for the WAG run. Run `codeml tmp0001.ctl` to obtain the updated `rst*` outputs.
3. Rename the appropriate codeml output `rst2` to `in.BV` in this directory.
4. Set `usedata = 2` in `mcmctree.ctl` and rerun `mcmctree mcmctree.ctl` for final dating. The included control file uses this setting. The resulting time tree is `FigTree.tre`.

## CAFE expansion and contraction

The full tree spans a time range that prevented effective convergence in this analysis. The dated tree and gene-family table were therefore split into Ascomycota (`asco`) and Basidiomycota (`basi`) and analyzed separately. Prepare these subtrees from the final dated tree before running CAFE; subtree extraction is not automated here.

Run the preparation scripts from `CAFE`:

1. **`1_family_file_filter.py`** maps the gene-count table's species columns through the shared trait file, adds the CAFE description column, removes `Total`, and excludes families with 100 or more genes in any species. Tables are saved in `cafe_related_files`.
2. Remove the NEXUS wrapper and dating annotations from `FigTree.tre`, retaining its dated branch lengths, and prepare `asco/tree_asco.txt` and `basi/tree_basi.txt`. **`2_scale_tree.py`** multiplies branch lengths by 100. Its default is the ascomycete subtree; change the input/output names to process the basidiomycete subtree as well.
3. **`3_split_genefamily.py`** selects the species columns matching each subtree and writes `asco/gene_family_asco.txt` and `basi/gene_family_basi.txt`.
4. Run CAFE separately in each subtree directory:

   ```sh
   # From CAFE/asco
   cafe5 -i gene_family_asco.txt -t tree_asco_scaled.txt -p -k 2 -o k2 -c 8

   # From CAFE/basi
   cafe5 -i gene_family_basi.txt -t tree_basi_scaled.txt -p -k 2 -o k2 -c 8
   ```

5. Run the postprocessing scripts from the corresponding `k2` directory. `1_clean_tree2fig.py` removes family counts and significance marks from a selected CAFE reconstruction saved as `tree2fig.nwk` (used to view the node numbering assigned by the software). `2_extract_significant_genes.sh` selects families significant at the family level and then extracts contractions at the nodes listed in `node.txt`. `3_extract_genes.sh` extracts changes without that significance filter: expansions in the supplied `asco` script and contractions in the `basi` script. These scripts copy the corresponding OrthoFinder family FASTA files from the configured sequence directory.

Positive values in `Gamma_change.tab` indicate expansions and negative values indicate contractions.

`cafe_related_files` contains the input tables and trees, while `timetree` and the two `k2` directories hold dating and CAFE run outputs.

## Software

The workflow uses OrthoFinder, Python with pandas, Biopython and ETE3, Perl with AnyDBM_File and Term::ANSIColor, PAML (`mcmctree` and `codeml`), and CAFE5. The optional taxonomy helper requires an internet connection and an Entrez email address.
