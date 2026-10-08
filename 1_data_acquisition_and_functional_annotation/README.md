# Data acquisition and functional annotation

These scripts summarize genome features, compare general yeasts with multicellular fungi, and annotate protein sequences. Only a few example inputs are included.

## Genome features

Run the numbered scripts in order from this directory, with the gene-family steps below completed before step 5. Adjust the input paths in the scripts for your working environment. Scripts that read species traits use `related_files/76_fungi_trait.txt`.

1. **`1_size_n_gc.py`** reads the genome FASTA files in `76_fna` and writes genome size (Mb) and GC content (%) to `genome_stat.tsv`.
2. **`2_te.sh`** builds a repeat library for each genome with RepeatModeler, then runs RepeatMasker. **`3_extract_te2tsv_masked.py`** reads the `bases masked` percentage from the RepeatMasker `.tbl` files in `76_fna` and adds it as `TE content(%)` to the summary table.
3. **`4_add_genecount.py`** uses the accession/name mappings to standardize species names, counts sequences in the corresponding protein FASTA files, and adds gene counts in thousands. It writes `genome_stat_updated.tsv`.
4. **`5_add_genefamily_count.py`** adds the number of orthogroups present in each species from `genefamily_count/Orthogroups.presence_counts.tsv`. The output is `genome_stat_updated_with_gf_counts.tsv`, with the column `Gene family number`.
5. **`6_plot_violin.py`** joins the genome summary with the species trait table and plots the five features. Each plot shows a two-sided Mann–Whitney U test comparing general yeasts (`trait = 1`) with multicellular fungi (`trait = 0`). These tests do not account for phylogeny.
6. **`7_pgls.py`** compares the same features while accounting for phylogenetic relatedness, as described below.

## Gene-family counts

Run these scripts from `genefamily_count` before adding gene-family counts to the genome summary:

1. **`1_create_count_table.py`** replaces the species column names in the OrthoFinder gene-count table using the species mapping and removes the `Total` column.
2. **`2_genefamily_loss_sig.py`** counts orthogroups with at least one gene in each species. A one-sided Mann–Whitney U test asks whether yeasts have fewer gene families. If significant, the script lists orthogroups absent from all yeasts but present in at least half of the multicellular species. This is a presence/absence screen, not an ancestral loss reconstruction.
3. **`3_plot_violin.py`** plots the per-species counts and annotates the one-sided comparison.

## Phylogenetic comparison

`7_pgls.py` reads the completed genome summary and rooted tree from `related_files`, and uses the `names` and `trait` columns in `related_files/76_fungi_trait.txt` for group assignments.

For each feature, the model is `feature ~ trait`, with Brownian-motion covariance and Pagel's lambda estimated by maximum likelihood between 0 and 1. The coefficient represents general yeast minus multicellular fungi. The script checks taxon matching and branch lengths, omits missing feature values separately for each model, and applies Benjamini–Hochberg correction across the five PGLS tests.

Run `python3 7_pgls.py`. Results, the merged data, a taxon-matching report, and PNG/PDF/SVG plots are saved in `related_files`. The plots show the observed distributions, annotated with PGLS p-values and adjusted q-values; stars reflect the q-values.

## Protein functional annotation

`parallel_annotate.sh` runs eggNOG-mapper in DIAMOND mode for the `.faa` files in `76_faa`, using the fungal database, and saves one set of annotation outputs per protein file in `76_anno`. Set the database location and CPU/job limits before running. This step can run independently of the genome-feature comparisons.

## Software

Python scripts use pandas, NumPy, Biopython, SciPy, matplotlib, and seaborn. Repeat analysis requires RepeatModeler and RepeatMasker; functional annotation requires eggNOG-mapper and its database. Both shell workflows use GNU Parallel.

## Example files

Example `genome_stat*.tsv` tables are in `related_files/`. Steps 1–6 use the current working directory; `7_pgls.py` reads from `related_files/`.
