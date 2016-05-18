# Code to calculate node rankings in the network
# It is assumed that the working directory is the source file location
############################################################################################################

############################################################################################################
#### Load libraries ####
# Clear R console screen output
cat("\014")

library(synapseClient)
library(knitr)
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
cl = makeCluster(14)
registerDoParallel(cl)

# Login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'rankNodesInModules.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", ref="branch", refName='rankNodes')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0(thisFileName))

# Synapse specific parameters
activityName = 'Node rankings in modules'
activityDescription = 'Ranking genes based on network structure'
############################################################################################################

############################################################################################################
#### Utility functions block ####

# Function to rank the nodes of network based on its enrichment for neighborhoods
scoreVertices <- function(g, order = 1){
  # g = unweighted, undirected graph with Gx as a vertex property
  # Gx = gene/vertex importance score
  # order = defines network neighborhood (i.e., include order number of nodes away from selected node)
  Gx = igraph::V(g)$Gx
  names(Gx) = igraph::V(g)$name
  
  neighbor.graph = igraph::make_ego_graph(g, order = order, nodes = igraph::V(g), mode = "all")
  names(neighbor.graph) = igraph::V(g)$name
  
  D = igraph::distances(g)
  
  S = igraph::strength(g, mode = "all")
  
  Nx = foreach(vert=V(g)$name,
               .combine = c,
               .packages = c('igraph'),
               .export = c('neighbor.graph', 'Gx', 'D', 'S'),
               .verbose = T) %dopar% {
                 
                 SP = igraph::shortest_paths(neighbor.graph[[vert]], from = vert, 
                                             to = V(neighbor.graph[[vert]]), mode = "all", 
                                             predecessors = T)
                 
                 Nx = Gx[igraph::V(neighbor.graph[[vert]])$name]
                 Nx = Nx*(1/D[vert,igraph::V(neighbor.graph[[vert]])$name])
                 Nx = Nx*(1/S[vert])
                 names(Nx) = V(neighbor.graph[[vert]])$name
                 sum(Nx[names(Nx) != vert])
               }
  names(Nx) = V(g)$name
  return(Nx)
}
############################################################################################################

############################################################################################################
#### Get differential expression and variants results from synapse ####
ALL_USED_IDs = c('syn5745241', 'syn5752526')

dexp = downloadFile('syn5745241') %>%
  dplyr::rename(ensembl_gene_id = V1) %>%
  dplyr::mutate(diffexp.score = -log10(adj.P.Val)*abs(logFC))
dexp$diffexp.score[dexp$diffexp.score < 0] = 0

variant = read.table(synGet('syn5752526')@filePath, sep = '\t', header = TRUE) %>%
  dplyr::select(Corrected.p.value, Gene.ensGene, AAChange.ensGene) %>%
  plyr::ddply(.(AAChange.ensGene), .fun = function(x){
    data.frame(Corrected.p.value = x$Corrected.p.value, 
               Gene.ensGene = unlist(str_split(x$Gene.ensGene, ',')))
  }) %>%
  dplyr::select(-AAChange.ensGene) %>%
  dplyr::rename(variant.adj.P.Val = Corrected.p.value, ensembl_gene_id = Gene.ensGene) %>%
  dplyr::mutate(variant.score = -log10(variant.adj.P.Val)) %>%
  unique
variant$variant.score[variant$variant.score < 0] = 0

Gx = full_join(dexp, variant) %>%
  unique
Gx$variant.adj.P.Val[is.na(Gx$variant.adj.P.Val)] = 1
Gx$variant.adj.P.Val[Gx$variant.adj.P.Val > 1] = 1
Gx$adj.P.Val[is.na(Gx$adj.P.Val)] = 1
Gx[is.na(Gx)] = 0

Gx = Gx %>%
  dplyr::mutate(gene.score = variant.score + diffexp.score) %>%
  unique() %>%
  group_by(ensembl_gene_id) %>%
  top_n(1, gene.score) %>%
  data.frame
rownames(Gx) = Gx$ensembl_gene_id
############################################################################################################

############################################################################################################
#### Get network ####
ALL_USED_IDs = c(ALL_USED_IDs, 'syn5907918')

# Get network from synapse (rda format)
load(synGet('syn5907918')@filePath)

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(bicNetworks$rankConsensus$network, mode = 'undirected', diag = F)

# Set Gx as vertex property
V(g)$Gx = Gx[V(g)$name,'gene.score']
V(g)$Gx[is.na(V(g)$Gx)] = 0
############################################################################################################

############################################################################################################
#### Get modules ####
# Download modules from synapse
modules = downloadFile('syn5923906')
ALL_USED_IDs = c(ALL_USED_IDs, 'syn5923906')

# Convert ensemble gene id's to hgnc symbols using biomart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = modules$EnsembleID, mart = ensembl)
############################################################################################################

############################################################################################################
#### Score/Rank genes ####

#### Based on neighborhood enrichment ####
vertex.score = scoreVertices(g, order = 3) %>% 
  rownameToFirstColumn('ensembl_gene_id') %>% 
  dplyr::rename(net.score = DF)
stopCluster(cl)

#### Based on node degree ####
node.degree = degree(g, mode = 'all') %>% 
  rownameToFirstColumn('ensembl_gene_id') %>% 
  dplyr::rename(degree.score = DF)
############################################################################################################

############################################################################################################
#### Store results in synapse ####
results = inner_join(vertex.score, node.degree) %>% 
  inner_join(Gx) %>%
  dplyr::select(ensembl_gene_id, adj.P.Val, logFC, variant.adj.P.Val, diffexp.score, variant.score, gene.score, net.score, degree.score) %>%
  left_join(modules %>%
              dplyr::rename(ensembl_gene_id = EnsembleID)) %>%
  left_join(ensg2hgnc) 

write.table(results, file = 'vertexRankings.tsv', row.names = F, sep = '\t', quote=F)
obj = File('vertexRankings.tsv', name = 'BIC Rank Consensus', parentId = 'syn6100411')
annotations(obj) = list(dataType = 'analysis', 
                        fileType = 'tsv',
                        normalizationStatus = 'TRUE',
                        analysisType = 'vertexRanking',
                        method = 'bic',
                        disease	= 'craniosynostosis',
                        organism = 'HomoSapiens')
obj = synStore(obj, activityName = activityName, 
               activityDescription = activityDescription, 
               used = ALL_USED_IDs, 
               executed = thisFile)

tmp = results %>%
  dplyr::mutate(ranks1 = rank(net.score),
                ranks2 = rank(rank(degree.score) * rank(net.score)),
                ranks3 = rank(rank(diffexp.score) * rank(variant.score) * rank(net.score)),
                ranks4 = rank(rank(diffexp.score) + rank(variant.score) + rank(net.score))) %>%
  dplyr::arrange(desc(ranks1))