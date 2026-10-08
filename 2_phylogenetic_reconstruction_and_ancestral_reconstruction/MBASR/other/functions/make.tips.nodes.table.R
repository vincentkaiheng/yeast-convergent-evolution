make.tips.nodes.table <- function(file.name.trait.data,file.name.tree,file.name.asr.results,multiplier) {

nodes.table=make.nodes.table(file.name.asr.results,multiplier)
tips.table=make.tips.table(file.name.trait.data,file.name.tree)

tips.nodes.table=rbind(tips.table,nodes.table)
tips.nodes.table=noquote(tips.nodes.table)

return(tips.nodes.table) }

