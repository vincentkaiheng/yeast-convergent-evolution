make.nodes.table <- function(file.name.asr.results,multiplier) {

asr.table=as.matrix(read.table(file.name.asr.results,sep="\t"))
rownames(asr.table)=asr.table[,1]
asr.table=asr.table[,-1]
mode(asr.table)="numeric"
colnames(asr.table)=NULL

n.states=dim(asr.table)[2]

if (n.states==2) { results=bin.nodes.from.reconstruction.2.states(asr.table,multiplier) }
if (n.states==3) { results=bin.nodes.from.reconstruction.3.states(asr.table,multiplier) }
if (n.states==4) { results=bin.nodes.from.reconstruction.4.states(asr.table,multiplier) }
if (n.states==5) { results=bin.nodes.from.reconstruction.5.states(asr.table,multiplier) }
if (n.states==6) { results=bin.nodes.from.reconstruction.6.states(asr.table,multiplier) }
if (n.states==7) { results=bin.nodes.from.reconstruction.7.states(asr.table,multiplier) }
if (n.states==8) { results=bin.nodes.from.reconstruction.8.states(asr.table,multiplier) }
if (n.states==9) { results=bin.nodes.from.reconstruction.9.states(asr.table,multiplier) }
if (n.states==10) { results=bin.nodes.from.reconstruction.10.states(asr.table,multiplier) }

return(results) }

