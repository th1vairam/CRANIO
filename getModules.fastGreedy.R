# Function to get bicNetworks from synapse and find modules and push results back to synapse

# Load libraries
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(igraph)
library(metanetwork)

library(synapseClient)
library(knitr)
library(githubr)

synapseLogin()

# Set run parameters
max.mod.sz = 500
min.mod.sz = 20

# Get latest commit of the executable from github (for provenance)
thisRepo <- githubr::getRepo(repository = "th1vairam/CRANIO", ref="branch", refName='moduleAnal')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath= 'getModules.fastGreedy.R')

# Get all bicNetworks.rda from the source project
bic.id = 'syn7342818'

# Get network from synapse
net.obj = synapseClient::synGet(bic.id)
load(net.obj@filePath)

# Get the adjacency matrix
adj <- bicNetworks$network

# Get modules
g = graph_from_adjacency_matrix(adj, mode = 'upper', weighted = NULL, diag = FALSE)
mod = igraph::fastgreedy.community(g)

# Get individual clusters from the igraph community object
geneModules = igraph::membership(mod) %>%
  unclass %>%
  as.data.frame %>%
  plyr::rename(c('.' = 'moduleNumber'))

geneModules = cbind(data.frame(Gene.ID = rownames(geneModules)),
                    geneModules)              

# Iteratively split modules that of maximum size mx.sz
count = 1
max.mod.num = max(geneModules$moduleNumber)
while(count <= max.mod.num){
  if (sum(geneModules$moduleNumber == count) > max.mod.sz){
    sg = induced_subgraph(g, V(g)$name %in% geneModules$Gene.ID[geneModules$moduleNumber == count])
    sg.mod = igraph::fastgreedy.community(sg) %>% membership
    geneModules[names(sg.mod), 'moduleNumber'] = sg.mod + max.mod.num
    max.mod.num = max(geneModules$moduleNumber)
  }
  count = count + 1
  cat(count);cat('\n')
}

# Rename modules with size less than min module size to 0
filteredModules = geneModules %>% 
  dplyr::group_by(moduleNumber) %>%
  dplyr::summarise(counts = length(unique(Gene.ID))) %>%
  dplyr::filter(counts >= min.mod.sz)
geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
geneModules$moduleNumber = as.factor(geneModules$moduleNumber) %>% unclass %>% as.numeric()
geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)

# Find modularity (Q) of the network
Q <- metanetwork::compute.Modularity(adj, geneModules, method = 'Newman1')
NQ <- metanetwork::compute.LocalModularity(adj, geneModules)
Qds <- metanetwork::compute.ModularityDensity(adj, geneModules)
# MQ <- metanetwork::compute.ModuleQualityMetric(adj, geneModules)

# Write results to synapse
algo = 'igraph.fast_greedy'
write.table(geneModules, file = paste0(algo,'.modules.tsv'), row.names=F, quote=F, sep = '\t')
obj = synapseClient::File(paste0(algo,'.modules.tsv'), name = paste('Modules',algo), parentId = net.obj@properties$parentId)
synapseClient::annotations(obj) = list('algorithm' = algo,
                                       'Q' = Q,
                                       'NQ' = NQ,
                                       'Qds' = Qds,
                                       'fileType' = 'tsv',
                                       'resultsType' = 'modules')
obj = synapseClient::synStore(obj, used = net.obj, executed = thisFile, activityName = 'Module Identification')