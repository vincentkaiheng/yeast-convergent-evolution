# Comparative genomic analyses of independent yeast lineages

This repository contains analysis scripts, example data and selected results for the study of independent yeast lineages.

## Analyses

| Directory | Contents |
|---|---|
| [1. Data and annotation](1_data_acquisition_and_functional_annotation/README.md) | Genome features, repeats and functional annotation |
| [2. Phylogeny](2_phylogenetic_reconstruction_and_ancestral_reconstruction/README.md) | Species trees and ancestral states |
| [3. Gene families](3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/README.md) | OrthoFinder, PAPS and CAFE |
| [4. Functional overlap](4_analysis_of_functional_convergence_between_shared_lost_genes_and_novel_genes/README.md) | Shared OGs, GO terms and PFAMs |
| [5. Relative rates](5_relative_evolutionary_rate_analysis/README.md) | RERconverge analyses |
| [6. Machine learning](6_machine_learning/README.md) | XGBRF, SHAP and rank aggregation |

## Running the scripts

Each folder has a README with the analysis steps and software requirements. Run scripts from the folder specified there and adjust the paths for your own setup. Example file locations are noted in the code, relative to the script folder.

Only part of the input data is included. Module 1 has a genome, proteome and annotation example for *Aureobasidium melanogenum* P16. Other folders contain example tables, trees and selected results.

The gene-family summaries in module 1 and the analyses in modules 4 and 6 use OrthoFinder results from module 3.

## Software sources

MBASR, PAPS and REVIGO sources are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). MrBayes should be installed separately. The MBASR license is included in its folder.
