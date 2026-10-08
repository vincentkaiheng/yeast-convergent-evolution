clear.the.mb.folder <- function() {

keep.these=NULL

all.files=dir()

if(length(all.files)>0) {

the.mb.exe=which(all.files=="mb.exe")
the.mb=which(all.files=="mb")

mb.exe.test=length(the.mb.exe)
mb.test=length(the.mb)

if(mb.exe.test==1 & mb.test==1) { keep.these=c(the.mb.exe,the.mb) }

if(mb.exe.test==1 & mb.test==0) { keep.these=c(the.mb.exe) }

if(mb.exe.test==0 & mb.test==1) { keep.these=c(the.mb) }

removal.list=all.files
if(length(keep.these)>0) { removal.list=removal.list[-keep.these] }

file.remove(removal.list)

}

if(length(all.files)==0) {
msg="FAILED. The \"mb\" folder is missing the MrBayes v3.2.7a executable file."
setwd(main.directory)
stop(msg)
}

if(length(keep.these)==0) {
msg="FAILED. The \"mb\" folder is missing the MrBayes v3.2.7a executable file."
setwd(main.directory)
stop(msg)
}

return(invisible(NULL)) }
