make.MrBayes.constraints.topology.and.node.ages <- function(file.name.tree,constraint.type) {

type_test=0
if(constraint.type=="hard") { type_test=type_test+1 }
if(constraint.type=="partial") { type_test=type_test+1 }
if(type_test==0) { stop("constraint.type must be \"hard\" or \"partial\"") }

suppressWarnings(suppressMessages(library(ape)))
library(ape)

suppressWarnings(suppressMessages(library(phytools)))
library(phytools)

tree=read.tree(file.name.tree)

n_internal_nodes=tree$Nnode
n_tips=length(tree$tip.label)
root_node_number=n_tips+1
last_node=root_node_number+n_internal_nodes-1

tip_numbers=seq(from=1,to=n_tips,by=1)
tip_names=tree$tip.label

internal_node_numbers=seq(from=root_node_number,to=last_node,by=1)
tree$node.label=internal_node_numbers

first_constraint=internal_node_numbers[1]
my_end=length(internal_node_numbers)
last_constraint=internal_node_numbers[my_end]
my_names=internal_node_numbers
node_names=paste("node",my_names,sep="")

results_list=list()
count_node=first_constraint
count_rep=1
repeat {
descendant_tips=getDescendants(tree,count_node)
remove_these=which(descendant_tips>n_tips)
if(length(remove_these)>0) { descendant_tips=descendant_tips[-remove_these] }
results_list[[count_rep]]=descendant_tips
count_node=count_node+1
count_rep=count_rep+1
if(count_node==last_constraint+1) break }
names(results_list)=node_names
#results_list

n_constraints=length(results_list)
constraint_results=list()
constraint_count=1
repeat {
current_node_name=node_names[constraint_count]
current_descendant_tip_numbers=results_list[[constraint_count]]
current_descendant_tip_names=tip_names[current_descendant_tip_numbers]
string_tip_names=paste(current_descendant_tip_names,collapse=" ")
x=paste("constraint",current_node_name,constraint.type,"=",string_tip_names,";",sep=" ")
constraint_results[[constraint_count]]=x
constraint_count=constraint_count+1
if(constraint_count==n_constraints+1) break }
constraint_results=unlist(constraint_results)
#constraint_results

head="prset topologypr = constraints ("
foot=");"
y=paste(node_names,collapse=",")
enforce_constraints=paste(head,y,foot,sep="")
#enforce_constraints

topological_constraints=c(constraint_results,enforce_constraints)
#topological_constraints

writeLines(topological_constraints,"temp1.txt",sep="\r\r")
temp1=readLines("temp1.txt",warn=F)
file.remove("temp1.txt")

last_internal_node_number=last_node
my_node_heights=nodeHeights(tree)
tallest_tip=max(my_node_heights)

results_1=list()
count_nodes=root_node_number
count=1

repeat {
node_of_interest=count_nodes
height_of_interest=nodeheight(tree,node_of_interest)
age_of_interest=height_of_interest-tallest_tip
age_of_interest=age_of_interest*-1
results_1[[count]]=age_of_interest
count_nodes=count_nodes+1
count=count+1
if(count_nodes==last_internal_node_number+1) break }

results_1=unlist(results_1)
results_1=round(results_1,digits=8)
root_age=results_1[1]
internal_node_ages=results_1[2:length(results_1)]

results_2=list()
count=1

repeat {
height_of_interest=nodeheight(tree,count)
age_of_interest=height_of_interest-tallest_tip
age_of_interest=age_of_interest*-1
results_2[[count]]=age_of_interest
count=count+1
if(count==n_tips+1) break }

results_2=unlist(results_2)
results_2=round(results_2,digits=8)
tip_ages=results_2

node_names=node_names[-1]

tip_age_constraints=paste(tip_names,"=fixed(",tip_ages,")",sep="")
node_age_constraints=paste(node_names,"=fixed(",internal_node_ages,")",sep="")
head2="calibrate"
foot2=";"
age_constraints=c(head2,tip_age_constraints,node_age_constraints,foot2)

x1="prset treeagepr=fixed("
x2=");"
x3=paste(x1,root_age,x2,sep="")

age_constraints=c(age_constraints,"",x3,"")

top_and_BLs=c(temp1,age_constraints)
writeLines(top_and_BLs,"MrBayes.topological.and.node.age.constraints.txt")

msg="MrBayes topological and node age constraints were written to file."
msg=noquote(msg)

return(msg) }

