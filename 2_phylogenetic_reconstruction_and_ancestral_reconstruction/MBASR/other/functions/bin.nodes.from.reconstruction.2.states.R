bin.nodes.from.reconstruction.2.states <- function(asr.table,multiplier) {

n.states=dim(asr.table)[2]
highest.state=n.states-1
n.nodes=dim(asr.table)[1]
state.seq=seq(from=0,to=highest.state,by=1)

results="0"

count=1
repeat {
current.probs=asr.table[count,]
names(current.probs)=state.seq
current.probs=sort(current.probs,decreasing=T)
winner.prob=current.probs[1]
runner.up.prob=current.probs[2]
mult.runner.up.prob=runner.up.prob*multiplier

if(winner.prob>=mult.runner.up.prob) { winner.state=names(winner.prob) }
if(winner.prob<mult.runner.up.prob) { winner.state="01" }

results=c(results,winner.state)

count=count+1
if(count==n.nodes+1) break }
results=results[-1]

node.names=rownames(asr.table)
node.states=results
results.table=cbind(node.names,node.states)
colnames(results.table)=NULL

results.table=noquote(results.table)

return(results.table) }

