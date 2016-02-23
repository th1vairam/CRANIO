source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")

require(synapseClient)
synapseLogin()

sparrowNetworkObj <- synGet('syn5652273',downloadLocation='./')
sparrowNetwork <- data.table::fread(sparrowNetworkObj@filePath,stringsAsFactors=FALSE,data.table=F)
rownames(sparrowNetwork) <- sparrowNetwork$V1
sparrowNetwork <- sparrowNetwork[,-1]



bar<-synGetActivity(sparrowNetworkObj)
source('grabCRANIOData.R')

baz <- read.csv('cranioRNAseq.csv',stringsAsFactors=F)
S <- cor(baz)
foo <- data.matrix(sparrowNetwork)[which(upper.tri(data.matrix(sparrowNetwork)))]
thresVal <- sort(foo^2,decreasing=T)[4e4]

sparrowNetwork <- sparrowNetwork^2
edgeList <- which(sparrowNetwork >=thresVal,T)
edval <- rep(0,nrow(edgeList))
for (i in 1:nrow(edgeList)){
  edval[i] <- sparrowNetwork[edgeList[i,1],edgeList[i,2]]
}
edgeList <- cbind(edgeList,edval)
colnames(edgeList) <- c('node1','node2','weight')
rownames(edgeList) <- paste0('e',1:nrow(edgeList))
edgeList <- data.frame(edgeList,stringsAsFactors=F)
library(dplyr)
edgeList <- arrange(edgeList,desc(weight))

source('covarianceSelectionBisection.R')
ggm <- covarianceSelectionBisection(S,rankedEdges=edgeList[,1:2],numberObservations=nrow(baz)*ncol(baz),lowerBoundEdge=1e4,upperBoundEdge=4e4)
