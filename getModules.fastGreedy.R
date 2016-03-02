#!/usr/bin/env Rscript

# Function to get modules from network adjacency matrix (from synapse as rda file)
# Get arguments from comman line
args = c('syn5700531')

# Clear R console screen output
cat("\014")

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)

# Needs the dev branch
library(rGithubClient)

# login to synapse
synapseLogin()

# Get github links for provenance
thisFileName = 'getModules.fastGreedy.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", 
                    ref="branch", 
                    refName='moduleAnal')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

# Synapse specific parameters
activityName = 'Module Identification'
activityDescription = 'Clustering network genes in to modules using fast greedy algorithm of igraph'

# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet(args[1])
FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name)
parentId = NET_OBJ$properties$parentId

# Load sparse network
load(NET_OBJ@filePath)

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(sparseNetwork, mode = 'undirected', weighted = T, diag = F)

# Get modules using fast.greedy method (http://arxiv.org/abs/cond-mat/0408187)
mod = igraph::fastgreedy.community(g)
collectGarbage()

# Get individual clusters from the igraph community object
clust.numLabels = igraph::membership(mod)
collectGarbage()

# Change cluster number to color labels
labels = WGCNA::labels2colors(clust.numLabels)

# Get results
geneModules = data.frame(GeneIDs = V(g)$name,
                         moduleNumber = as.numeric(clust.numLabels), 
                         modulelabels = labels)

# Write results to synapse
algo = str_replace_all(algorithm(mod), '[^[:alnum:]]', '_')

write.table(geneModules, paste(FNAME,algo,'tsv',sep='.'), sep='\t', row.names=F, quote=F)
MOD_OBJ = File(paste(FNAME,algo,'tsv',sep='.'), name = paste(FNAME,algo,'Modules'), parentId = parentId)
annotations(MOD_OBJ) = annotations(NET_OBJ)
MOD_OBJ@annotations$fileType = 'tsv'
MOD_OBJ@annotations$moduleMethod = paste('igraph',algo,sep=':')
MOD_OBJ@annotations$modularity = modularity(mod)
MOD_OBJ = synStore(MOD_OBJ, 
                   executed = thisFile,
                   used = NET_OBJ,
                   activityName = activityName,
                   activityDescription = activityDescription)