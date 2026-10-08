plot.tree.with.pie.charts <- function(file.name.tree,file.name.trait.data,file.name.plot.settings) {

file.name.asr="MrBayes.ASR.results.txt"

source(file.name.plot.settings)

state.decision.multiplier=state.decision.odds

my.tip.label.offset=my.tip.label.shift.right
my.tip.pie.offset=my.tip.pie.shift.right

suppressWarnings(suppressMessages(library(ape)))
library(ape)

suppressWarnings(suppressMessages(library(phytools)))
library(phytools)

tree=read.tree(file.name.tree)

n_tips=length(tree$tip.label)
n_nodes=tree$Nnode
root_node_number=n_tips+1
last_node=root_node_number+n_nodes-1
node_numbers=seq(from=root_node_number,to=last_node,by=1)

tree$node.label=node_numbers

BLs.orig=tree$edge.length
n.edges=length(BLs.orig)
mean.BL=mean(BLs.orig)

scaling.factor=1
if(mean.BL<0.1) { scaling.factor=100 }
if(mean.BL<0.01) { scaling.factor=1000 }
if(mean.BL<0.001) { scaling.factor=10000 }
if(mean.BL<0.0001) { scaling.factor=100000 }
if(mean.BL<0.00001) { scaling.factor=1000000 }
if(mean.BL<0.000001) { scaling.factor=10000000 }

BLs.scaled=BLs.orig*scaling.factor
tree$edge.length[1:n.edges]=BLs.scaled

ASR_table=read.table(file.name.asr,sep="\t")
ASR_table=as.matrix(ASR_table)
rownames(ASR_table)=ASR_table[,1]
ASR_table=ASR_table[,-1]
mode(ASR_table)="numeric"
colnames(ASR_table)=NULL

tip_data=read.table(file.name.trait.data,sep="\t")
tip_data=as.matrix(tip_data)
colnames(tip_data)=NULL
rownames(tip_data)=tip_data[,1]
tip_data=as.matrix(tip_data[,-1])
colnames(tip_data)="state"

xyz=strsplit(tip_data[,1],"")
xyz=unlist(xyz)
options(warn=-1)
mode(xyz)="numeric"
options(warn=0)
max_state=max(xyz,na.rm=T)

first_state=0
all_states=seq(from=first_state,to=max_state,by=1)
n_states=length(all_states)

needed_rows=dim(tip_data)[1]
needed_cols=length(all_states)

state_names=seq(from=0,to=max_state,by=1)
state_names=paste("state ",state_names,sep="")

c0="red"
c1="orange"
c2="yellow"
c3="green"
c4="cyan"
c5="blue"
c6="purple"
c7="magenta"
c8="gray60"
c9="gray90"
my.colors=c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)

my.colors=my.colors[1:n_states]

tips_matrix=matrix(0,nrow=needed_rows,ncol=needed_cols)
rownames(tips_matrix)=rownames(tip_data)

current_taxon=1
repeat {
current_score=tip_data[current_taxon,1]
test_and=grep("&",current_score)

if(length(test_and)>0) {
current_score=strsplit(current_score,"&")
current_score=unlist(current_score)
options(warn=-1)
mode(current_score)="numeric"
options(warn=0)
destination_columns=current_score+1
current_strength=1/length(current_score)
round(current_strength,digits=3)
tips_matrix[current_taxon,destination_columns]=current_strength
}

if(length(test_and)==0) {
options(warn=-1)
mode(current_score)="numeric"
options(warn=0)
destination_columns=current_score+1
current_strength=1
tips_matrix[current_taxon,destination_columns]=current_strength
}

current_taxon=current_taxon+1
if(current_taxon==needed_rows+1) break }

Qs=which(tip_data=="?")
if(length(Qs)>0) {
tips_matrix[Qs,1:needed_cols]="?"
}

#tips_matrix

tips_list_correct_order=tree$tip.label
tips_list_incorrect_order=rownames(tips_matrix)
my_match=match(tips_list_incorrect_order,tips_list_correct_order)

temp_matrix=cbind(tips_matrix,tips_list_correct_order)
temp_dim=dim(temp_matrix)
temp_cols=temp_dim[2]

abc=colnames(temp_matrix)
abc[temp_cols]="correct_seq"
colnames(temp_matrix)=abc

temp_matrix[,temp_cols]=my_match
options(warn=-1)
mode(temp_matrix)="numeric"
options(warn=0)

temp_matrix=as.matrix(temp_matrix[order(temp_matrix[,temp_cols]),])

dim_again=dim(tips_matrix)
cols_again=dim_again[2]
rows_again=dim_again[1]

tips_matrix[1:rows_again,1:cols_again]=temp_matrix[1:rows_again,1:cols_again]
options(warn=-1)
mode(tips_matrix)="numeric"
options(warn=0)

pdf(file="tree.plot.pdf",width=my.plot.width,height=my.plot.height)

if(show.node.numbers==FALSE) {
plot.phylo(tree,direction="rightwards",no.margin=TRUE,cex=my.tip.label.font.size,label.offset=my.tip.label.offset,edge.width=my.tree.edge.width,underscore=T)
}

if(show.node.numbers==TRUE) {
plot.phylo(tree,direction="rightwards",no.margin=TRUE,cex=my.tip.label.font.size,label.offset=my.tip.label.offset,edge.width=my.tree.edge.width,underscore=T,show.node.label=T)
}

par(lwd=my.node.pie.border.width)
nodelabels(pie=ASR_table,cex=my.node.pie.size,piecol=my.colors)

par(lwd=my.tip.pie.border.width)
tiplabels(pie=tips_matrix,cex=my.tip.pie.size,adj=my.tip.pie.offset,piecol=my.colors)

legend(x=0,y=my.legend.shift.up,legend=state_names,box.col="white",title="",cex=my.legend.size,fill=my.colors,horiz=T)

if(show.time.scale==TRUE) {

all.dist.nodes=dist.nodes(tree)
root.dist.nodes=all.dist.nodes[,root_node_number]
root.tips.dists=root.dist.nodes[1:n_tips]
root.age=max(root.tips.dists)

my.time.scale.n.ticks=10

one.tick.size=root.age/(my.time.scale.n.ticks-1)
scale.seq=seq(from=0,to=root.age,by=one.tick.size)

scale.seq.2=pretty(scale.seq,n=my.time.scale.n.ticks)
new.seq.length=length(scale.seq.2)
scale.seq.2[new.seq.length]=root.age

tick.labels=rev(scale.seq.2)
short.int=tick.labels[1]-tick.labels[2]
norm.int=tick.labels[2]-tick.labels[3]
adj.val=norm.int-short.int

scale.seq.3=scale.seq.2-adj.val
scale.seq.3[1]=scale.seq.2[1]
scale.seq.3[new.seq.length]=scale.seq.2[new.seq.length]

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xscale <- range(lastPP$xx)
tscale <- c(0,root.age)
beta <- diff(xscale)/diff(tscale)
alpha <- xscale[1] - beta * tscale[1]
ticks.at=beta*scale.seq.3+alpha

tick.labels[1]=NA

tick.labels=tick.labels+my.time.scale.youngest.tip
#tick.labels=round(tick.labels,digits=my.time.scale.label.n.digits)

tick.labels=tick.labels/scaling.factor

starting.options=options()
starting.scipen=starting.options$scipen
options(scipen=999)

axis(side=1,pos=my.time.scale.bar.shift.up,at=ticks.at,labels=tick.labels,cex.axis=my.tip.label.font.size,lwd=my.tree.edge.width,padj=-1*my.time.scale.label.shift.up,tcl=-0.1)

options(scipen=starting.scipen)

}

if(show.edge.transitions==TRUE) {

multiplier=state.decision.multiplier
file.name.asr.results="MrBayes.ASR.results.txt"
tips.nodes.table=make.tips.nodes.table(file.name.trait.data,file.name.tree,file.name.asr.results,multiplier)
edge.transition.table=make.edge.transition.table(file.name.tree,tips.nodes.table)
edge.transition.table[,4]=edge.transition.table[,4]+my.time.scale.youngest.tip
edge.transition.table[,5]=edge.transition.table[,5]+my.time.scale.youngest.tip
my.edge.labels=edge.transition.table[,8]

my.edge.transitions.shift.up=my.edge.transitions.shift.up*-1
edgelabels(my.edge.labels,frame="none",cex=my.edge.transitions.font.size,adj=c(my.edge.transitions.shift.left,my.edge.transitions.shift.up),col="#888888")
}

dev.off()

msg="Tree plot was written to file."
msg=noquote(msg)

return(msg) }

