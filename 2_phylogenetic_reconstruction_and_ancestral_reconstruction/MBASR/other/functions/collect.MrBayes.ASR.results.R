collect.MrBayes.ASR.results <- function(file.name.p) {

report.sd="no"
#report.sd="yes"

my_options=options()
my_starting_scipen=my_options$scipen
my_starting_digits=my_options$digits
options(scipen=999)
options(digits=12)

temp_data=readLines(file.name.p)
temp_data=temp_data[-1]
writeLines(temp_data,"temp.txt")

data=read.table("temp.txt",header=T)
file.remove("temp.txt")
last_row=dim(data)[1]

xx=as.matrix(data[,1])
yy=as.matrix(data[,2])

LnLs=as.matrix(data[,2])
remove.these.A=seq(from=1,to=501,by=1)
LnLs.2=as.matrix(LnLs[-remove.these.A,])
mean.LnL=mean(LnLs.2)
my_threshold=mean.LnL

max.LnL=max(LnLs)
min.LnL=min(LnLs)

xxx=c(xx[1],xx[length(xx)])
yyy=c(my_threshold,my_threshold)

xxxx=c(50000,50000)
yyyy=c(min.LnL,max.LnL)

pdf("gens_v_likelihood.pdf",width=10,height=5.625)
plot(xx,yy,type="n",xlab="Generation",ylab="LnL")
points(xx,yy,pch=20,cex=0.25)
points(xxx,yyy,type="l",col="red")
points(xxxx,yyyy,type="l",col="red")
dev.off()

data=data[-remove.these.A,]
my_keep=which(data[,2]>my_threshold)
keep_gens=data[my_keep,]

n.keep=length(my_keep)
write.table(n.keep,"n.keep.txt",quote=F,col.names=F,row.names=F)

first_needed_col=8
n_cols=dim(keep_gens)[2]
keep_gens=keep_gens[,first_needed_col:n_cols]
n_cols=dim(keep_gens)[2]

my_names=names(keep_gens)
my_names=gsub("\\.\\.1\\.","_",my_names)
my_names=gsub("\\.","",my_names)
names(keep_gens)=my_names

my_means=colMeans(keep_gens)

#standard deviation calculations
all_sds=0
sd_count=1
repeat {
sd_data=keep_gens[,sd_count]
my_sd=sd(sd_data,na.rm=T)
all_sds=c(all_sds,my_sd)
sd_count=sd_count+1
if(sd_count==n_cols+1) break }
all_sds=all_sds[-1]
my_sds=all_sds

test_highest_state=my_names[1:11]
test_highest_state=strsplit(test_highest_state,"_node")
test_highest_state=matrix(unlist(test_highest_state),ncol=2,byrow=T)
test_highest_state=test_highest_state[,1]
test_highest_state=gsub("p","",test_highest_state)
mode(test_highest_state)="numeric"
test_highest_state=max(test_highest_state)
#test_highest_state
needed_n_cols=test_highest_state+1

state_names=seq(from=0,to=test_highest_state,by=1)
state_names=paste("state_",state_names,sep="")

node_names=strsplit(my_names,"_")
node_names=matrix(unlist(node_names),ncol=2,byrow=T)
node_names=node_names[,2]
qwerty=seq(from=needed_n_cols,to=length(node_names),by=needed_n_cols)
node_names=node_names[qwerty]

x=matrix(my_means,ncol=needed_n_cols,byrow=T)
colnames(x)=state_names
rownames(x)=node_names
mode(x)="numeric"
x=round(x,digits=8)

#make sd matrix
sd_matrix=matrix(my_sds,ncol=needed_n_cols,byrow=T)
colnames(sd_matrix)=state_names
rownames(sd_matrix)=node_names
mode(sd_matrix)="numeric"
sd_matrix=round(sd_matrix,digits=8)

my_end=dim(x)[1]
count=1
repeat {
my_sum=sum(x[count,])
my_dif=1-my_sum
my_max=max(x[count,])
my_max_pos=which(x[count,]==my_max)
x[count,my_max_pos]=my_max+my_dif
count=count+1
if(count==my_end+1) break }

write.table(format(x,digits=8),"MrBayes.ASR.results.txt",col.names=F,row.names=T,sep="\t",quote=F)

if(report.sd=="yes") {
write.table(format(sd_matrix,digits=8),"standard.deviations.txt",col.names=F,row.names=T,sep="\t",quote=F)
file.copy("standard.deviations.txt",results.directory)
file.remove("standard.deviations.txt")
}

options(scipen=my_starting_scipen)
options(digits=my_starting_digits)

msg="MrBayes ASR results were wrtitten to file."
msg=noquote(msg)

return(msg) }

