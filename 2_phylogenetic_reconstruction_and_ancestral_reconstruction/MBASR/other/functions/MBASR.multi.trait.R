MBASR.multi.trait <- function(file.name.tree,file.name.multi.trait.csv,file.name.plot.settings,file.name.ordered.characters,n.samples) {

setwd(input.files.directory)
convert.csv.matrix.to.single.traits(file.name.multi.trait.csv)

setwd(multi.traits.directory)
all.files=dir()
n.files=length(all.files)

setwd(input.files.directory)

test=readLines(file.name.ordered.characters,warn=F)
mode(test)="numeric"
options(warn=-1)
test=max(test,na.rm=T)
options(warn=0)
if(test>0) { test=1 }
if(test=="-Inf") { test=0 }

if(test>0) {
options(warn=-1)
ordered.characters=read.table(file.name.ordered.characters)
options(warn=0)
ordered.characters=unlist(ordered.characters)
names(ordered.characters)=NULL
mode(ordered.characters)="numeric"
}

character.types=rep("unordered",n.files)
if(test>0) { character.types[ordered.characters]="ordered" }

msg=""
write.table(msg,quote=F,row.names=F,col.names=F)

count=1
repeat {

msg=paste("Trait: ",count,sep="")
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

current.character.type=character.types[count]

setwd(multi.traits.directory)
current.file.name.trait.data=all.files[count]
file.copy(current.file.name.trait.data,input.files.directory)

setwd(main.directory)
MBASR(file.name.tree,current.file.name.trait.data,file.name.plot.settings,current.character.type,n.samples)

setwd(input.files.directory)
file.remove(current.file.name.trait.data)

temp.name=gsub(".txt","",current.file.name.trait.data)
table.name=paste(temp.name,"_ASR.results.txt",sep="")
pdf.name=paste(temp.name,"_tree.plot.pdf",sep="")
setwd(results.directory)
file.rename("MrBayes.ASR.results.txt",table.name)
file.rename("tree.plot.pdf",pdf.name)

setwd(main.directory)

count=count+1
if(count==n.files+1) break }

setwd(multi.traits.directory)
all.files.again=dir()
new.file.names=gsub(".txt","_tree.plot.pdf",all.files.again)
first.file.name=new.file.names[1]
last.file.name=new.file.names[length(new.file.names)]
last.file.number=gsub("trait_","",last.file.name)
last.file.number=gsub("_tree.plot.pdf","",last.file.number)
first.file.number=gsub("trait_","",first.file.name)
first.file.number=gsub("_tree.plot.pdf","",first.file.number)
combined.name=paste("_traits_",first.file.number,"-",last.file.number,"_tree.plots.pdf",sep="")

setwd(results.directory)
pdf_combine(new.file.names,output=combined.name)

setwd(input.files.directory)
unlink("multi.traits",recursive=T)
setwd(main.directory)

return(invisible(NULL)) }

