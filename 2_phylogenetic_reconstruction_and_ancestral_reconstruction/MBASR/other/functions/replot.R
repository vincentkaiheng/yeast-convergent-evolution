replot <- function(file.name.tree,file.name.trait.data,file.name.plot.settings) {

#main.directory=MBASR.directory
#functions.directory=paste(main.directory,"/functions",sep="")
#input.files.directory=paste(main.directory,"/input.files",sep="")
#mb.directory=paste(main.directory,"/mb",sep="")
#results.directory=paste(main.directory,"/results",sep="")

setwd(results.directory)
file.copy("MrBayes.ASR.results.txt",input.files.directory)

setwd(input.files.directory)

plot.tree.with.pie.charts(file.name.tree,file.name.trait.data,file.name.plot.settings)

file.remove("MrBayes.ASR.results.txt")
file.copy("tree.plot.pdf",results.directory,overwrite=T)
file.remove("tree.plot.pdf")

setwd(main.directory)

msg1="Tree replot was written to file."
msgX=""
msg=c(msgX,msg1,msgX)
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

return(invisible(NULL)) }
