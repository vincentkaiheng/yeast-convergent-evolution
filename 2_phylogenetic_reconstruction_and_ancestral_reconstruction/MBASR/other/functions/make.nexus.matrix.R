make.nexus.matrix <- function(my_matrix){

temp=my_matrix
temp=apply(temp,1,paste,collapse="")
temp=as.matrix(temp)

my_n_taxa=length(temp)

all_scores=temp[,1]
score_chars=strsplit(all_scores,"")
score_lengths=lengths(score_chars)
my_n_chars=min(score_lengths)

tax=paste("NTAX=",my_n_taxa,sep="")
char=paste("NCHAR=",my_n_chars,sep="")

line1="#NEXUS"
line2=""
line3="BEGIN DATA;"
line4=paste("DIMENSIONS ",tax," ",char,";",sep="")
line5="FORMAT DATATYPE = STANDARD GAP = - MISSING = ?;"
line6="MATRIX"
line7=""

line_n1=";"
line_n2="END;"

taxa_names=rownames(temp)
ntaxa=length(taxa_names)

count=1
current_taxon=taxa_names[[count]]
current_scores=temp[[count]]
blank_line=""
current_lines=rbind(current_taxon,current_scores,blank_line)
rownames(current_lines)=NULL
growing_lines=current_lines

count=2
repeat {
current_taxon=taxa_names[[count]]
current_scores=temp[[count]]
blank_line=""
current_lines=rbind(current_taxon,current_scores,blank_line)
rownames(current_lines)=NULL
growing_lines=rbind(growing_lines,current_lines)
count=count+1
if(count==ntaxa+1) break }

header=rbind(line1,line2,line3,line4,line5,line6,line7)
footer=rbind(line_n1,line_n2)
rownames(header)=NULL
rownames(footer)=NULL

my_matrix=rbind(header,growing_lines,footer)

write.table(my_matrix,"matrix.nex",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,na="")

msg="Nexus matrix was written to file."
msg=noquote(msg)

return(msg) }

