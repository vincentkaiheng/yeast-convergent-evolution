make.edge.transition.table <- function(file.name.tree,tips.nodes.table) {

tree=read.tree(file.name.tree)

edge.matrix=tree$edge
n.edges=dim(edge.matrix)[1]
my.seq=seq(from=1,to=n.edges,by=1)
edge.matrix=cbind(edge.matrix,edge.matrix)
edge.matrix=edge.matrix[,-1]
edge.matrix[,1]=my.seq
my.col.names=c("edge.number","start.node","end.node")
colnames(edge.matrix)=my.col.names

edge.matrix=cbind(edge.matrix,edge.matrix,edge.matrix)
edge.matrix=edge.matrix[,-9]
edge.matrix[,4:8]=NA
my.col.names=c("edge.number","start.node","end.node","start.age","end.age","start.state","end.state","state.transition")
colnames(edge.matrix)=my.col.names

dist.mat=dist.nodes(tree)
n.tips=length(tree$tip.label)
root.node=n.tips+1
root.to.tips.dists=dist.mat[,root.node]
max.height=max(root.to.tips.dists)
node.ages=root.to.tips.dists-max.height
node.ages=node.ages*-1
node.ages=round(node.ages,digits=6)
node.ages=as.matrix(node.ages)

count=1
repeat{
current.node=edge.matrix[count,2]
current.node.age=node.ages[current.node,1]
edge.matrix[count,4]=current.node.age
count=count+1
if(count==n.edges+1) break }

count=1
repeat{
current.node=edge.matrix[count,3]
current.node.age=node.ages[current.node,1]
edge.matrix[count,5]=current.node.age
count=count+1
if(count==n.edges+1) break }

edge.matrix=as.data.frame(edge.matrix)

count=1
repeat{
current.node=edge.matrix[count,2]
current.node.state=tips.nodes.table[current.node,2]
edge.matrix[count,6]=current.node.state
count=count+1
if(count==n.edges+1) break }

count=1
repeat{
current.node=edge.matrix[count,3]
current.node.state=tips.nodes.table[current.node,2]
edge.matrix[count,7]=current.node.state
count=count+1
if(count==n.edges+1) break }

edge.matrix[,8]=""

count=1
repeat{
current.start.state=edge.matrix[count,6]
current.end.state=edge.matrix[count,7]
identical.test=identical(current.start.state,current.end.state)
if(identical.test==FALSE) {
 #tag=paste(count,": ",current.start.state," -> ",current.end.state,sep="")
 tag=paste(current.start.state," -> ",current.end.state,sep="")
 edge.matrix[count,8]=tag
 }
count=count+1
if(count==n.edges+1) break }

blank.these=grep("\\?",edge.matrix[,8])
if(length(blank.these)>0) { edge.matrix[blank.these,8]="" }

return(edge.matrix) }

