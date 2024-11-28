## function ##

pres_abs <- function(metadata, roary_output,last_column,output){

metadata <- metadata
roary_output <- roary_output
last_column <- as.numeric(last_column)
output <- output
  
m <- read.csv(metadata, header=TRUE, sep="\t")
dim(m)
head(m)

dat <- read.csv(roary_output, header=TRUE, check.names = FALSE)
dim(dat)
length(names(dat))
length(unique(names(dat)))
names(dat)[1:30]

head(dat[1:10,1:25])
dbt <- dat[,c(1,15:ncol(dat))]
dim(dbt)
x <- gsub("[[:alnum:]].*","1",as.matrix(dbt[,2:ncol(dbt)]))
y <- gsub("$.*","0",x)
z <- gsub("10","1",y)
w <- data.frame(Gene=dbt$Gene,z)
v <- as.data.frame(t(w))
p <- v[c(2:nrow(v)),]
names(p) <- v[1,]
class(p)
dim(p)
row(p)

q <- data.frame(sample=names(dbt[,2:ncol(dbt)]),p)
dim(q)
q[1:30,1:10]

head(m)
names(m) <- c("sample",names(m)[2:last_column])
q4 <- merge(m,q,by="sample",all.x=FALSE)
dim(q4)
names(q4)[1:20]

print(q4[1:10,1:25])

q4 <- return(q4)

write.table(q4,output, row.names=FALSE, quot=F)

}
