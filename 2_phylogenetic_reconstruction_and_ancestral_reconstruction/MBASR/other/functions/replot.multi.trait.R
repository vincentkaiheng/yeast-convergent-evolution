replot.multi.trait <- function(file.name.tree,file.name.multi.trait.csv,file.name.plot.settings) {

setwd(input.files.directory)
convert.csv.matrix.to.single.traits(file.name.multi.trait.csv)

setwd(multi.traits.directory)
all.files=dir()
n.files=length(all.files)

setwd(multi.traits.directory)

count=1
repeat {
current.trait.data.file.name=all.files[count]
file.copy(current.trait.data.file.name,input.files.directory)
current.results.file.name=gsub(".txt","_ASR.results.txt",current.trait.data.file.name)
setwd(results.directory)
file.copy(current.results.file.name,input.files.directory)
setwd(input.files.directory)
file.rename(current.results.file.name,"MrBayes.ASR.results.txt")
plot.tree.with.pie.charts(file.name.tree,current.trait.data.file.name,file.name.plot.settings)
file.remove("MrBayes.ASR.results.txt")
file.remove(current.trait.data.file.name)
file.copy("tree.plot.pdf",results.directory)
file.remove("tree.plot.pdf")
setwd(results.directory)
current.pdf.file.name=gsub(".txt","_tree.plot.pdf",current.trait.data.file.name)
file.rename("tree.plot.pdf",current.pdf.file.name)
setwd(multi.traits.directory)
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

msgX=""
msgZ="Tree replots were written to file."
write.table(msgX,row.names=F,col.names=F,quote=F)
flush.console()
write.table(msgZ,row.names=F,col.names=F,quote=F)
flush.console()
write.table(msgX,row.names=F,col.names=F,quote=F)
flush.console()

return(invisible(NULL)) }

