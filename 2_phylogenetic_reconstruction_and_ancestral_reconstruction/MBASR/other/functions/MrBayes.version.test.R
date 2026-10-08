MrBayes.version.test <- function() {

correct.mb.version=FALSE

OS.test=Sys.info()
OS.test=grep("Windows",OS.test)
OS.test=length(OS.test)

if(OS.test==1) { my_executable="./mb.exe" }
if(OS.test==0) { my_executable="./mb" }

options(warn=-1)
system2(my_executable,"mb.version.test.txt",stdout=tempfile("stdout.txt"))
options(warn=0)

version.log=readLines("version.test.txt")

my.test=grep("3.2.7",version.log)
my.test.2=length(my.test)

if(my.test.2==1) { correct.mb.version=TRUE }

return(correct.mb.version) }

