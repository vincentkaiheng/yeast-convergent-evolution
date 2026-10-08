repair.root.multi.singles <- function(file.name.tree) {

t1=read.tree(file.name.tree)

binary.test=is.binary(t1)
if(binary.test==FALSE) {

t2=t1
file.name.tree.orig=paste(file.name.tree,"_original.nwk",sep="")
file.rename(file.name.tree,file.name.tree.orig)

t2=multi2di(t2,random=F)
t2=collapse.singles(t2,root.edge=TRUE)
write.tree(t2,file.name.tree)
}

return(invisible(NULL)) }

