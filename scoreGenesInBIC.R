# Score and rank nodes of CRANIO (case+control) BIC network 
## It is assumed that the working directory is the source file location
############################################################################################################

############################################################################################################
#### Load libraries ####
# Clear R console screen output
cat("\014")

library(synapseClient)
library(githubr)

library(CovariateAnalysis) # Refere th1vairam repo in github
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(igraph)
library(biomaRt)

library(doParallel)
library(foreach)

cl = makeCluster(6)
registerDoParallel(cl)

# Login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'scoreGenesInBIC.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", ref="branch", refName='rankNodes')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0(thisFileName))

# Synapse specific parameters
activityName = 'Score nodes'
activityDescription = 'Score nodes based on network structure'
############################################################################################################

############################################################################################################
#### Utility functions block ####
# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genes.in.significant.set, # A character vector of differentially expressed or some significant genes to test
                             genes.in.gene.set, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genes.in.background # Background genes 
){
  genes.in.significant.set = intersect(genes.in.significant.set, genes.in.background) # back ground filtering
  genes.in.nonsignificant.set = base::setdiff(genes.in.background, genes.in.significant.set)
  genes.in.gene.set = intersect(genes.in.gene.set, genes.in.background) # back ground filtering
  genes.our.gene.set = base::setdiff(genes.in.background,genes.in.gene.set)
  
  tp = length(intersect(genes.in.gene.set, genes.in.significant.set))             
  fp = length(intersect(genes.in.gene.set, genes.in.nonsignificant.set))
  fn = length(intersect(genes.our.gene.set, genes.in.significant.set))
  tn = length(intersect(genes.our.gene.set, genes.in.nonsignificant.set))
  
  pval = fisher.test(matrix(c(tp, fp, fn, tn),nrow=2, ncol=2), alternative="greater")
  odds = (tp*tn)/(fp*fn)
  
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genes.in.gene.set),
                    noverlap = tp,
                    odds = odds,
                    Genes = paste(intersect(genes.in.gene.set, genes.in.significant.set), collapse = '|')
  )
  )
}

# Function to calculate node scores
score.nodes <- function(g, G, h = 3, mode = 'all'){
  
  # Get background genes from graph
  background.genes = igraph::V(g)$name
  G = G[intersect(background.genes, names(G))]
  G.not = rep(0, length(setdiff(background.genes, names(G))))
  names(G.not) = setdiff(background.genes, names(G))
  G = c(G, G.not)
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = mode, mindist = hi)
  }, g)
  
  # Find summed weight of current node based on h-layer neighbors
  node.scores = foreach::foreach (x = 1:length(neighbor.nodes), .combine = cbind) %dopar% {
    score = sapply(neighbor.nodes[[x]], function(y, G){
      sum(G[names(y)] * (1/node.degree[names(y)]))
    }, G)
  } 
  node.scores = rowSums(node.scores * t(matrix(rep(1/c(1:h), length(G)), h, length(G))), na.rm = T)
  node.scores = node.scores + G
  names(node.scores) = background.genes
  
  return(node.scores)
}

# Function to find weighted regulators on directed networks
regulatorAnalysis.undirected_weighted <- function(g, G, h = 3, n = 100, FDR = 0.05){
  
  # Convert adjacency to igraph object
  background.genes = igraph::V(g)$name
  
  # Compute node scores based on neighborhood scores
  node.scores = score.nodes(g, G, h, mode = 'all')
  
  # Permute node labels and calculate null scores
  perm.node.scores = foreach::foreach(i = 1:n, 
                                      .combine = cbind, 
                                      .packages = c('igraph', 'dplyr', 'parallel', 'doParallel', 'foreach'),
                                      .export = c('score.nodes'),
                                      .verbose = TRUE) %dopar% {
                                        pg = igraph::permute(g, sample(1:length(background.genes), length(background.genes)))
                                        perm.node.scores = score.nodes(pg, G, h, mode = 'all')
                                      } 
  
  # Perform one sample t-test to estimate significance
  pval = foreach::foreach(i = 1:length(background.genes), .combine = rbind) %dopar% {
    tmp = t.test(perm.node.scores[i,], mu = node.scores[i], alternative = 'less')
    data.frame(pval = tmp$p.value, t = tmp$statistic, t.low = tmp$conf.int[1], t.high = tmp$conf.int[2])
  }
  fdr = p.adjust(pval$pval, method = 'fdr')
  names(fdr) = background.genes
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find leaf nodes to identify global regulators
  node.in.degree = igraph::degree(g, mode = 'all')
  
  # Promote high degree nodes as global regulators 
  key.regulators = list()
  key.regulators$scores = node.scores
  key.regulators$fdr = fdr
  key.regulators$regulators = background.genes[(fdr <= FDR)]
  key.regulators$global.regulators = background.genes[((node.degree > (mean.node.degree + 2*stddev.node.degree)) | (node.in.degree == 0)) & (fdr <= FDR)]
  key.regulators$local.regulators = setdiff(key.regulators$regulators, key.regulators$global.regulators)
  
  return(key.regulators)
}
############################################################################################################

############################################################################################################
#### Get cranio case+control BIC network from synapse ####
net.id = 'syn7342818'
all.used.ids = net.id
load(synGet(net.id)@filePath)

g = igraph::graph_from_adjacency_matrix(bicNetworks$network, mode = 'upper', diag = F)
############################################################################################################

############################################################################################################
#### Get differential expression and variants results from synapse ####
diff.exp.ids = c(case.control = 'syn5745241')
all.used.ids = c(all.used.ids, diff.exp.ids)
de.genesets = llply(diff.exp.ids, function(de.id){
  de = CovariateAnalysis::downloadFile(de.id) %>%
    dplyr::mutate(Score = abs(logFC) * (-1) * log10(adj.P.Val)) %>%
    daply(.(V1), .fun = function(x){x$Score})
})

variants.id = 'syn5752526'
all.used.ids = c(all.used.ids, variants.id)
variants = downloadFile(variants.id) %>%
  dplyr::select(one_of('Gene.ensGene', 'Corrected p-value')) %>%
  plyr::rename(c('Gene.ensGene' = 'Gene.ID', 'Corrected p-value' = 'adj.P.Value')) %>%
  daply(.(Gene.ID), function(x){
    sum(-log10(x$adj.P.Value))
  })
  
all.genesets = c(de.genesets, list(variants = variants))
############################################################################################################

############################################################################################################
#### Score genes in network ####
all.genesets = lapply(all.genesets, function(x,g){
  x = x[intersect(names(x), V(g)$name)]
},g)

all.scores = lapply(all.genesets, function(x, g){
  tmp = regulatorAnalysis.undirected_weighted(g, x, h = 3, n = 5, FDR = 0.05)
  
  tmp1 = join_all(list(tmp$scores %>%
                         CovariateAnalysis::rownameToFirstColumn('Gene.ID') %>%
                         dplyr::rename(Scores = DF),
                       tmp$fdr %>%
                         CovariateAnalysis::rownameToFirstColumn('Gene.ID') %>%
                         dplyr::rename(FDR = DF)),
                  type = 'full')
  
  tmp1$regulator = FALSE
  tmp1$regulator[tmp1$Gene.ID %in% tmp$regulator] = TRUE
  
  tmp1$global.regulator = FALSE
  tmp1$global.regulator[tmp1$Gene.ID %in% tmp$global.regulator] = TRUE
  
  tmp1$local.regulator = FALSE
  tmp1$local.regulator[tmp1$Gene.ID %in% tmp$local.regulator] = TRUE
  
  return(tmp1)
}, g)
############################################################################################################

############################################################################################################
#### Store results in synapse ####
objs = lapply(names(all.scores), function(listName, allScores){
  write.table(allScores[[listName]], file = paste0('GeneScores',listName,'.tsv'), row.names = F, sep = '\t', quote=F)
  obj = File(paste0('GeneScores',listName,'.tsv'), 
             name = paste('Gene Scores',listName), 
             parentId = 'syn7320952')
  annotations(obj) = list(fileType = 'tsv',
                          resultsType = 'geneScores',
                          algorithm = 'weighted.undirected.regulators',
                          dataName = 'cranio')
  obj = synStore(obj, activityName = activityName, 
                 activityDescription = activityDescription, 
                 used = as.character(all.used.ids), 
                 executed = thisFile)
}, all.scores)

stopCluster(cl)
############################################################################################################