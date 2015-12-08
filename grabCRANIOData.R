###

require(synapseClient)
synapseLogin()

foo <- synGet('syn4586436')
require(data.table)

bar <- fread(foo@filePath,data.table=F,stringsAsFactors = F)
n1 <- scan(foo@filePath,what='character',nlines=1)
dim(bar)

n1 <- c('GeneID',n1)
colnames(bar) <- n1
bar[1:5,1:5]

require(dplyr)
bar[1:5,1:5]
hist(bar[,6])

cranioExpr <- bar[,-1]
cranioExpr <- apply(cranioExpr,2,as.numeric)
cranioExpr[1:5,1:5]
rownames(cranioExpr) <- bar[,1]
cranioExpr[1:5,1:5]
cranioExpr <- t(cranioExpr)

winsorize <- function(x,per=.99){
  up <- quantile(x,per)
  low <- quantile(x,1-per)
  x[x>=up] <- up
  x[x<=low] <- low
  return(x)
}
cranioExpr <- apply(cranioExpr,2,winsorize)
cranioExpr <- scale(cranioExpr)


write.csv(cranioExpr,file='cranioRNAseq.csv',quote=F)


