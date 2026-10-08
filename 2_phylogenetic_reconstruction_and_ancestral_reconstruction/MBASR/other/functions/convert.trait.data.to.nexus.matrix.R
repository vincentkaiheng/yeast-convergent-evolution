convert.trait.data.to.nexus.matrix <- function(file.name.trait.data) {

trait_data=read.table(file.name.trait.data,header=F,sep="\t")
trait_data=as.matrix(trait_data)
colnames(trait_data)=NULL

n_taxa=dim(trait_data)[1]
test_and=grep("&",trait_data)
if(length(test_and)>0) {
temp_cells=trait_data[test_and]
temp_cells=gsub("&"," ",temp_cells)
temp_cells=paste("(",temp_cells,")",sep="")
trait_data[test_and]=temp_cells
}

my_rownames=trait_data[,1]
my_scores=trait_data[,2]

my_matrix=matrix(,ncol=5,nrow=length(my_rownames))
rownames(my_matrix)=my_rownames

my_matrix[,1]=my_scores
my_matrix[,2]="?"
my_matrix[,3]="?"
my_matrix[,4]="?"
my_matrix[,5]="?"

make.nexus.matrix(my_matrix)

msg="Nexus matrix was written to file."
msg=noquote(msg)

return(msg) }

