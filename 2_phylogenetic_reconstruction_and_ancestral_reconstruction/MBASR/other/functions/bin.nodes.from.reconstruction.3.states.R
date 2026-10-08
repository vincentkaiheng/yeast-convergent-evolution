bin.nodes.from.reconstruction.3.states <- function(asr.table,multiplier) {

n.states=dim(asr.table)[2]
highest.state=n.states-1
n.nodes=dim(asr.table)[1]
state.seq=seq(from=0,to=highest.state,by=1)

results="0"
count=1

repeat {

done=0
current.probs=asr.table[count,]

names(current.probs)=state.seq
current.probs=sort(current.probs,decreasing=T)
winner.prob=current.probs[1]
runner.up.prob=current.probs[2]
mult.runner.up.prob=runner.up.prob*multiplier

 if(winner.prob>=mult.runner.up.prob) { 
 winner.state=names(winner.prob) 
 done=1
 }

 if(done==0) { 
 winner.prob=current.probs[2]
 runner.up.prob=current.probs[3]
 mult.runner.up.prob=runner.up.prob*multiplier

 if(winner.prob>=mult.runner.up.prob) {
 temp.states=names(current.probs)
 temp.states=temp.states[1:2]
 temp.states=sort(temp.states,decreasing=F)
 temp.states=paste(temp.states,sep="",collapse="")
 winner.state=temp.states
 done=1
 }
 }

 if(done==0) {
 winner.state="012"
 done=1
 }

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

