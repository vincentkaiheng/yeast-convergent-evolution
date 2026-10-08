convert.csv.matrix.to.single.traits <- function(file.name.csv) {

seq1=seq(from=1,to=9,by=1)
seq2=seq(from=10,to=99,by=1)
seq3=seq(from=100,to=999,by=1)
seq4=seq(from=1000,to=9999,by=1)

seq1b=paste("trait_","000",seq1,".txt",sep="")
seq2b=paste("trait_","00",seq2,".txt",sep="")
seq3b=paste("trait_","0",seq3,".txt",sep="")
seq4b=paste("trait_",seq4,".txt",sep="")

all.file.names=c(seq1b,seq2b,seq3b,seq4b)

data=read.csv(file.name.csv,header=F)
rownames(data)=data[,1]
data=data[,-1]
data=as.matrix(data)
colnames(data)=NULL
data=noquote(data)

temp=rownames(data)
temp=gsub(" ","_",temp)
rownames(data)=temp

n.chars=dim(data)[2]
if(n.chars>9999) { data=data[,1:9999] }

files.in.directory=dir()
multi.traits.folder.test=grep("multi.traits",files.in.directory)
multi.traits.folder.test=length(multi.traits.folder.test)
if(multi.traits.folder.test>0) { unlink("multi.traits",recursive=T) }

dir.create("multi.traits")
setwd(multi.traits.directory)

count=1
repeat {
current.col=as.matrix(data[,count])
outfile=all.file.names[count]
write.table(current.col,outfile,sep="\t",col.names=F,row.names=T,quote=F)
count=count+1
if(count==n.chars+1) break }

setwd(input.files.directory)

msg="All characters were indivdually written to file."
#write.table(msg,row.names=F,col.names=F,quote=F)
#flush.console()

return(invisible(NULL)) }

