validate.tree.and.trait.data <- function (file.name.tree,file.name.trait.data) {

test_tree=read.tree(file.name.tree)
test_data=read.table(file.name.trait.data,header=F,sep="\t")

tree_taxa=test_tree$tip.label
tree_taxa=as.matrix(tree_taxa)
tree_taxa=as.matrix(tree_taxa[order(tree_taxa[,1]),])

data_taxa=test_data[,1]
data_taxa=as.matrix(data_taxa)
data_taxa=as.matrix(data_taxa[order(data_taxa[,1]),])

n.tree.taxa=length(tree_taxa)
n.data.taxa=length(data_taxa)
culprit_taxa=c(setdiff(tree_taxa,data_taxa),setdiff(data_taxa,tree_taxa))
test_identical=identical(tree_taxa,data_taxa)
test_rooted=is.rooted(test_tree)

if(n.tree.taxa>n.data.taxa) {

msg="FAILED: There are MORE taxa in the tree than in the trait data. Please make corrections and try again."

setwd(main.directory)
stop(msg)
}

if(n.tree.taxa<n.data.taxa) {

msg="FAILED: There are LESS taxa in the tree than in the trait data. Please make corrections and try again."

setwd(main.directory)
stop(msg)
}

if(test_identical==FALSE) {

msg="FAILED: Taxon names in tree and trait data files are NOT identical. Please make corrections and try again."

setwd(main.directory)
stop(msg)
}

#if(test_rooted==FALSE) {

#msg="FAILED: Tree is NOT rooted. Please root the tree and try again."

#setwd(main.directory)
#stop(msg)
#}

tree.BLs=test_tree$edge.length
BL.test=length(tree.BLs)

if(BL.test==0) {
msg="FAILED: Tree has no branch lengths."

setwd(main.directory)
stop(msg)
}

test_data_2=test_data
rownames(test_data_2)=test_data_2[,1]
test_data_2[,1]=test_data_2[,2]
colnames(test_data_2)=NULL

character.string=test_data_2[,1]
character.string.1=gsub("&","",character.string)
character.string.2=gsub("\\?","",character.string.1)
character.string.3=gsub("0","",character.string.2)
character.string.4=gsub("1","",character.string.3)
character.string.5=gsub("2","",character.string.4)
character.string.6=gsub("3","",character.string.5)
character.string.7=gsub("4","",character.string.6)
character.string.8=gsub("5","",character.string.7)
character.string.9=gsub("6","",character.string.8)
character.string.10=gsub("7","",character.string.9)
character.string.11=gsub("8","",character.string.10)
character.string.12=gsub("9","",character.string.11)

character.string.13=paste(character.string.12,sep="",collapse="")
character.test=nchar(character.string.13)

if(character.test>0) {

character.string.14=gsub(""," ",character.string.13)

msg1="Legal trait characters are: 0 1 2 3 4 5 6 7 8 9 ? &"
msg2=paste("Illegal character(s) found:",character.string.14,sep="")
msg3="FAILED: There are illegal characters in your trait data. Please make corrections and try again."

write.table(msg1,quote=F,row.names=F,col.names=F)
flush.console()
write.table(msg2,quote=F,row.names=F,col.names=F)
flush.console()

setwd(main.directory)
stop(msg3)
}

return(invisible(NULL)) }

