MBASR <- function(file.name.tree,file.name.trait.data,file.name.plot.settings,character.type,n.samples) {

cleanup_mb_files="yes"
#cleanup_mb_files="no"

n.samples=validate.n.samples(n.samples)
n.samples.orig=n.samples
n.samples=n.samples+500

setwd(input.files.directory)
validate.tree.and.trait.data(file.name.tree,file.name.trait.data)

setwd(input.files.directory)
repair.root.multi.singles(file.name.tree)

setwd(input.files.directory)
prep.tree(file.name.tree)

setwd(mb.directory)
clear.the.mb.folder()

setwd(functions.directory)
file.copy("mb.prop.test.txt",mb.directory)
file.copy("mb.version.test.txt",mb.directory)
file.copy("mb.settings.txt",mb.directory)

setwd(mb.directory)
correct.mb.version=MrBayes.version.test()
if(correct.mb.version==FALSE) {
setwd(main.directory)
msg="FAILED: Incorrect MrBayes version number. Please use version 3.2.7 or 3.2.7a."
stop(msg)
}

my.starting.options=options()
my.starting.scipen=my.starting.options$scipen
options(scipen=999)

my_ngen=n.samples*100
my_ngen=paste("mcmcp ngen=",my_ngen,";",sep="")
my_ctype=character.type
my_ctype=paste("ctype ",my_ctype,": 1 ;",sep="")

options(scipen=my.starting.scipen)

setwd(mb.directory)
writeLines(my_ngen,"my_ngen.txt")
writeLines(my_ctype,"my_ctype.txt")

setwd(input.files.directory)

constraint.type="hard"
make.MrBayes.constraints.topology.and.node.ages("MBASR.prepped.tree.nwk",constraint.type)
file.copy("MrBayes.topological.and.node.age.constraints.txt",mb.directory)
file.remove("MrBayes.topological.and.node.age.constraints.txt")
file.remove("MBASR.prepped.tree.nwk")

setwd(input.files.directory)

convert.trait.data.to.nexus.matrix(file.name.trait.data)
file.copy("matrix.nex",mb.directory)
file.remove("matrix.nex")

OS.test=Sys.info()
OS.test=grep("Windows",OS.test)
OS.test=length(OS.test)

if(OS.test==1) { my_executable="./mb.exe" }
if(OS.test==0) { my_executable="./mb" }

msg="MrBayes is working on the ASR..."
msg=c("",msg)
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

setwd(mb.directory)

options(warn=-1)
system2(my_executable,"mb.prop.test.txt",stdout=tempfile("stdout.txt"))
options(warn=0)

make.proposal.exclusions("prop.test.txt")

options(warn=-1)
system2(my_executable,"mb.settings.txt",stdout=tempfile("stdout.txt"))
options(warn=0)

if(cleanup_mb_files=="yes") {
file.remove("my_ngen.txt")
file.remove("my_ctype.txt")
file.remove("my_props.txt")
file.remove("matrix.nex.t")
file.remove("matrix.nex.mcmc")
file.remove("matrix.nex")
file.remove("MrBayes.topological.and.node.age.constraints.txt")
file.remove("log.txt")

file.remove("mb.prop.test.txt")
file.remove("mb.settings.txt")
file.remove("mb.version.test.txt")

file.remove("prop.test.txt")
file.remove("version.test.txt")
}

msg="Summarizing the results..."
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

setwd(mb.directory)

file.name.p="matrix.nex.p"

collect.MrBayes.ASR.results(file.name.p)
n.keep=read.table("n.keep.txt")
n.keep=as.numeric(n.keep)

if(cleanup_mb_files=="yes") {
file.remove("matrix.nex.p")
file.remove("gens_v_likelihood.pdf")
file.remove("n.keep.txt")
}

file.copy("MrBayes.ASR.results.txt",input.files.directory)
file.remove("MrBayes.ASR.results.txt")

msg="Plotting the tree..."
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

setwd(input.files.directory)

plot.tree.with.pie.charts(file.name.tree,file.name.trait.data,file.name.plot.settings)

file.copy("MrBayes.ASR.results.txt",results.directory,overwrite=T)
file.copy("tree.plot.pdf",results.directory,overwrite=T)

file.remove("MrBayes.ASR.results.txt")
file.remove("tree.plot.pdf")

setwd(main.directory)

msg="Results were written to file."
msg=c(msg,"")
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

n.samples=n.samples.orig

n.discard=n.samples-n.keep

percent.discard=n.discard/n.samples*100
percent.discard=round(percent.discard,digits=1)
percent.summarized=100-percent.discard
percent.summarized=round(percent.summarized,digits=1)

sample.msg.1=paste("Samples requested: ",n.samples,sep="")
sample.msg.2=paste("Percent discarded: ",percent.discard,sep="")
sample.msg.3=paste("Percent summarized: ",percent.summarized,sep="")

write.table(sample.msg.1,row.names=F,col.names=F,quote=F)
flush.console()
write.table(sample.msg.2,row.names=F,col.names=F,quote=F)
flush.console()
write.table(sample.msg.3,row.names=F,col.names=F,quote=F)
flush.console()
msg=""
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

return(invisible(NULL)) }

