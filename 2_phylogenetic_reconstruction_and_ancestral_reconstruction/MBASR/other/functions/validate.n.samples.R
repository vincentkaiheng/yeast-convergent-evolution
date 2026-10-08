validate.n.samples <- function(n.samples) {

n.samples2=round(n.samples,digits=0)

if(n.samples2>10000) {
msg="FAILED: n.samples must be an integer ranging from 250 to 10000."
setwd(main.directory)
stop(msg)
}

if(n.samples2<250) {
msg="FAILED: n.samples must be an integer ranging from 250 to 10000. Please make the correction and try again."
setwd(main.directory)
stop(msg)
}

return(n.samples2) }

