setwd(MBASR.directory)

main.directory=MBASR.directory
functions.directory=paste(main.directory,"/other/functions",sep="")
input.files.directory=paste(main.directory,"/input",sep="")
mb.directory=paste(main.directory,"/other/mb",sep="")
results.directory=paste(main.directory,"/output",sep="")
multi.traits.directory=paste(input.files.directory,"/multi.traits",sep="")

msg=""
write.table(msg,quote=F,row.names=F,col.names=F)

msg="Toolkit: MBASR -- MrBayes Ancestral States with R (v2023.10.27)"
write.table(msg,quote=F,row.names=F,col.names=F)
flush.console()

#msg=""
#write.table(msg,quote=F,row.names=F,col.names=F)
#flush.console()

required.library.1="ape"
required.library.2="phytools"
required.library.3="pdftools"
cran.mirror="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"

my.libraries=installed.packages()
my.libraries=rownames(my.libraries)

test.1=which(my.libraries==required.library.1)
test.1=length(test.1)

test.2=which(my.libraries==required.library.2)
test.2=length(test.2)

test.3=which(my.libraries==required.library.3)
test.3=length(test.3)

all.required.libraries=sum(test.1,test.2,test.3)

if(all.required.libraries<3) {
msg="At least one dependency library is not yet installed."
write.table(msg,quote=F,row.names=F,col.names=F)
}

if(test.1==0) {
msg=paste("Installing package: ",required.library.1,sep="")
write.table(msg,quote=F,row.names=F,col.names=F)
install.packages(required.library.1,repos=cran.mirror,quiet=T)
}

if(test.2==0) {
msg=paste("Installing package: ",required.library.2,sep="")
write.table(msg,quote=F,row.names=F,col.names=F)
install.packages(required.library.2,repos=cran.mirror,quiet=T)
}

if(test.3==0) {
msg=paste("Installing package: ",required.library.3,sep="")
write.table(msg,quote=F,row.names=F,col.names=F)
install.packages(required.library.3,repos=cran.mirror,quiet=T)
}

msg=""
write.table(msg,quote=F,row.names=F,col.names=F)

msg="Loading dependency package: ape"
write.table(msg,quote=F,row.names=F,col.names=F)

suppressWarnings(suppressMessages(library(ape)))
library(ape)

msg="Loading dependency package: phytools"
write.table(msg,quote=F,row.names=F,col.names=F)

suppressWarnings(suppressMessages(library(phytools)))
library(phytools)

msg="Loading dependency package: pdftools"
write.table(msg,quote=F,row.names=F,col.names=F)

suppressWarnings(suppressMessages(library(pdftools)))
library(pdftools)

setwd(functions.directory)

all.files=dir()
exclude.1=which(all.files=="mb.prop.test.txt")
exclude.2=which(all.files=="mb.settings.txt")
exclude.3=which(all.files=="mb.version.test.txt")
exclude.these=c(exclude.1,exclude.2,exclude.3)
all.functions=all.files[-exclude.these]
n.functions=length(all.functions)
count=1
repeat {
source(all.functions[count])
count=count+1
if(count==n.functions+1) break }

setwd(main.directory)

msg=""
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

msg="Ready."
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

msg=""
write.table(msg,row.names=F,col.names=F,quote=F)
flush.console()

remove(all.files)
remove(exclude.1)
remove(exclude.2)
remove(exclude.these)
remove(all.functions)
remove(n.functions)
remove(count)

OS.test=Sys.info()
OS.test=grep("Windows",OS.test)
OS.test=length(OS.test)

if(OS.test==1) { my_exec="mb.exe" }
if(OS.test==0) { my_exec="mb" }

setwd(mb.directory)
files.in.mb.dir=dir()
setwd(main.directory)

executable.test=which(files.in.mb.dir==my_exec)
executable.test=length(executable.test)

if(executable.test==0) {
msg="WARNING: Please place an OS-specific executable of MrBayes 3.2.7a in the \"other/mb\" folder and rename it \"mb.exe\" (Windows) or \"mb\" (MacOS)." 
write.table(msg,quote=F,row.names=F,col.names=F)

msg=""
write.table(msg,quote=F,row.names=F,col.names=F)
}

remove(msg)
setwd(main.directory)

