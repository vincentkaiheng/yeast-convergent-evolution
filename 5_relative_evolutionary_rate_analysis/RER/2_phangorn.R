# Estimate branch lengths on the supplied species-tree topology.
library(RERconverge)

alignment_dir <- "shared_orthologs_sp_only"
tree_file <- "related_files/fungi_tree.tre"
output_file <- "related_files/gene_trees.tre"

estimatePhangornTreeAll(
    alndir = alignment_dir,
    treefile = tree_file,
    output.file = output_file
)
