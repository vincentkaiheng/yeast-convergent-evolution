make.proposal.exclusions <- function(log.file.name) {

my.log=readLines(log.file.name)

move.line.numbers=grep("-- Move",my.log)
move.lines=my.log[move.line.numbers]

move.lines.2=gsub(" ","",move.lines)

moves.lines.3=strsplit(move.lines.2,"=")
moves.lines.3=matrix(unlist(moves.lines.3),ncol=2,byrow=T)
moves.lines.3=moves.lines.3[,2]

multiplier.line.numbers=grep("Multiplier",moves.lines.3)
moves.lines.4=moves.lines.3[-multiplier.line.numbers]

n.moves=length(moves.lines.4)
moves.seq=seq(from=1,to=n.moves,by=1)

tau.lines.numbers=grep("Tau",moves.lines.4)

if(length(tau.lines.numbers)>0) { non.tau.lines.numbers=moves.seq[-tau.lines.numbers] }
if(length(tau.lines.numbers)==0) { non.tau.lines.numbers=moves.seq }

if(length(tau.lines.numbers)>0) { tau.lines=moves.lines.4[tau.lines.numbers] }
non.tau.lines=moves.lines.4[non.tau.lines.numbers]

if(length(tau.lines.numbers)>0) {
tau.props=strsplit(tau.lines,"\\(")
tau.props.2=matrix(unlist(tau.props),ncol=2,byrow=T)
tau.props.2=tau.props.2[,1]
}

non.tau.props=strsplit(non.tau.lines,"\\(")
non.tau.props.2=matrix(unlist(non.tau.props),ncol=2,byrow=T)
non.tau.props.2=non.tau.props.2[,1]

if(length(tau.lines.numbers)>0) { tau.post.text="(Tau{all},V{all})$prob=0;" }
non.tau.post.text="(V{all})$prob=0;"

if(length(tau.lines.numbers)>0) {
tau.final=paste("propset ",tau.props.2,tau.post.text,sep="")
non.tau.final=paste("propset ",non.tau.props.2,non.tau.post.text,sep="")
all.final=c(tau.final,non.tau.final)
}

if(length(tau.lines.numbers)==0) {
non.tau.final=paste("propset ",non.tau.props.2,non.tau.post.text,sep="")
all.final=non.tau.final
}

out.file.name="my_props.txt"
writeLines(all.final,out.file.name)

return(invisible(NULL)) }
