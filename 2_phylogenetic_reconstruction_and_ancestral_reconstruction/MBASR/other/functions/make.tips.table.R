make.tips.table <- function(file.name.trait.data,file.name.tree) {

traits=as.matrix(read.table(file.name.trait.data,sep="\t"))
rownames(traits)=traits[,1]
traits=traits[,-1]
colnames(traits)=NULL
traits=as.matrix(traits)
traits=noquote(traits)
traits=traits[sort(row.names(traits)),]
traits=as.matrix(traits)

n.taxa=dim(traits)[1]

count=1
repeat {
temp=traits[count,1]
temp=gsub("&","",temp)
n.chars=nchar(temp)
names(n.chars)=NULL
if(n.chars>1) { temp=strsplit(temp,"") }
if(n.chars>1) { temp=unlist(temp) }
if(n.chars>1) { temp=sort(temp) }
if(n.chars>1) { temp=paste(temp,collapse="") }
if(n.chars>1) { traits[count,1]=temp }
count=count+1
if(count==n.taxa+1) break }

tree=read.tree(file.name.tree)
tip.names=tree$tip.label
n.tips=length(tip.names)
tip.seq=seq(from=1,to=n.tips,by=1)

tip.matrix=cbind(tip.seq,tip.names)
rownames(tip.matrix)=tip.matrix[,2]
colnames(tip.matrix)=NULL
tip.matrix=tip.matrix[,-2]
tip.matrix=as.matrix(tip.matrix)
tip.matrix=noquote(tip.matrix)
tip.matrix=tip.matrix[sort(row.names(tip.matrix)),]
tip.matrix=as.matrix(tip.matrix)
tips=tip.matrix

tips.traits=cbind(tips,traits)
tips.traits=as.data.frame(tips.traits)
colnames(tips.traits)=NULL
tips.traits[,1]=as.numeric(tips.traits[,1])
tips.traits=tips.traits[order(tips.traits[,1]),]

terminal.names=paste("tip",tips.traits[,1],sep="")
terminal.states=tips.traits[,2]
terminal.matrix=cbind(terminal.names,terminal.states)
colnames(terminal.matrix)=NULL
results=noquote(terminal.matrix)

return(results) }

