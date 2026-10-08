prep.tree <- function(file.name.tree) {

scale.value=100

minBL=1/256

out.name="MBASR.prepped.tree.nwk"

starting.options=options()
starting.scipen=starting.options$scipen
options(scipen=999)

tree1=read.tree(file.name.tree)

tree2=tree1
tree2$node.label=NULL

node.heights=nodeHeights(tree2)
tree.height=max(node.heights)

scale.proportion=scale.value/tree.height

BLs.orig=tree2$edge.length

BLs.new=BLs.orig*scale.proportion

short.edges=which(BLs.new<minBL)

if(length(short.edges)>0) { BLs.new[short.edges]=minBL }

tree2$edge.length=BLs.new

write.tree(tree2,out.name)

options(scipen=starting.scipen)

return(invisible(NULL)) }

