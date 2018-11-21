# Code to calculate node rankings in the network
# It is assumed that the working directory is the source file location
############################################################################################################

############################################################################################################
#### Load libraries ####
# Clear R console screen output
cat("\014")

library(synapser)
library(knitr)
library(githubr)

library(CovariateAnalysis) # Refere th1vairam repo in github
library(data.table)
library(plyr)
library(tidyverse)

library(igraph)
library(biomaRt)
library(xlsx)

library(future)
library(furrr)

plan(multiprocess)

# Login to synapse
synLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'rankNodesInModulesTFDNA.R'

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
scoreVertices <- function(g, order = 2){
  # g = unweighted, undirected graph with Gx as a vertex property
  # Gx = gene/vertex importance score
  # order = defines network neighborhood (i.e., include order number of nodes away from selected node)
  Gx = igraph::V(g)$Gx
  names(Gx) = igraph::V(g)$name
  
  neighbor.graph = igraph::make_ego_graph(g, order = order, nodes = igraph::V(g), mode = "all")
  names(neighbor.graph) = igraph::V(g)$name
  
  D = igraph::distances(g)
  
  S = igraph::strength(g, mode = "all")
  
  Nx = map(V(g)$name,.f = function(vert){
    
    SP = igraph::shortest_paths(neighbor.graph[[vert]], from = vert, 
                                to = V(neighbor.graph[[vert]]), mode = "all", 
                                predecessors = T)
    
    Nx = Gx[igraph::V(neighbor.graph[[vert]])$name]
    Nx = Nx*(1/D[vert,igraph::V(neighbor.graph[[vert]])$name])
    Nx = Nx*(1/S[vert])
    names(Nx) = V(neighbor.graph[[vert]])$name
    sum(Nx[names(Nx) != vert])
  }) %>% 
    unlist()
  names(Nx) = V(g)$name
  return(Nx)
}
############################################################################################################

############################################################################################################
#### Get differential expression and variants results from synapse ####
downloadFile <- function(id){
  fread(synGet(id)$path, header = T, data.table = F)
}

all.used.ids = c('syn8555303', 'syn8555329')
dexp = c(cohort1 = 'syn8555303', cohort2 = 'syn8555329') %>%
  map(.f = function(id){
    downloadFile(id) %>%
      plyr::dlply(.(Model, Comparison), .fun = function(y){ 
        y %>% 
          dplyr::mutate(dscore = -log10(adj.P.Val) * abs(logFC)) %>%
          dplyr::select(ensembl_gene_id, hgnc_symbol, logFC, adj.P.Val, dscore)
      })
  })

## Get variants data from synapse
all.used.ids = c(all.used.ids, c('syn17015623', 'syn17015632'))
variants = c(all = 'syn17015623', any = 'syn17015632') %>%
  map(.f = function(id){
    downloadFile(id) %>%
      dplyr::select(Gene.refGene, `GERP++_RS`, CADD_phred, Polyphen2_HDIV_score, Polyphen2_HVAR_score) %>%
      dplyr::rename(GeneName = Gene.refGene,
                    GERP = `GERP++_RS`,
                    CADD = CADD_phred,
                    Polyphen2HDIV = Polyphen2_HDIV_score,
                    Polyphen2HVAR = Polyphen2_HVAR_score) %>%
      dplyr::mutate_at(c('GERP', 'CADD', 'Polyphen2HDIV', 'Polyphen2HVAR'), .funs = as.numeric) %>%
      dplyr::group_by(GeneName) %>%
      dplyr::mutate(rk = rank(rank(GERP) + rank(CADD) + rank(Polyphen2HVAR) + rank(Polyphen2HDIV))) %>%
      dplyr::top_n(1, rk) %>%
      dplyr::mutate(rk = rank(rank(GERP) + rank(CADD))) %>%
      dplyr::top_n(1, rk) %>%
      dplyr::mutate(rk = rank(rank(GERP))) %>%
      dplyr::top_n(1, rk) %>%
      distinct() %>%
      dplyr::select(-(rk))
  })

## Load modules from synapse
all.used.ids = c(all.used.ids, c('syn17015634', 'syn17015635'))
mod = c(cohort1 = 'syn17015634', cohort2 = 'syn17015635') %>%
  map(.f =function(tblId){
    synapser::synTableQuery(paste0('select name,id from ', tblId))$asDataFrame() %>%
      dplyr::filter(name == 'output.RData') %>%
      dplyr::select(id) %>%
      dplyr::mutate(folderName = id) %>%
      mutate_at(.vars = 'folderName', .funs = function(id){
        map(id, .f = function(x){
          synGet(x, downloadFile = F)$properties$parentId
        }) %>% 
          map(.f = function(x){
            synGet(x, downloadFile = F)$properties$name
          }) %>% 
          unlist()
      }) %>%
      group_by(folderName) %>%
      nest() %>%
      deframe() %>%
      map(.f = function(id){
        load(synGet(id)$path)
        
        # Remove small or big modules
        mod = output$modules
        sz = sapply(mod, length)
        mod = mod[sz >= 30 & sz <= 2000]
        return(mod)
      })
  })

## Important modules
imp.mod = list(
  cohort1 = list(control = c(),
                 coronal = c('c1_28', 'c1_3'),
                 metopic = c('c1_7', 'c1_54', 'c1_310', 'c1_363', 'c1_205', 'c1_59'),
                 sagittal = c('c1_47')),
  cohort2 = list(control = c('c1_13', 'c1_11', 'c1_32', 'c1_242', 'c1_123', 'c1_57', 
                             'c1_128', 'c1_114', 'c1_4', 'c1_177', 'c1_55'),
                 coronal = c(),
                 metopic = c('c1_101', 'c1_117', 'c1_146', 'c1_241', 'c1_152', 'c1_158'),
                 sagittal = c('c1_64')))

## Subset important modules
mod = mod %>%
  map2(imp.mod, .f = function(x,y){
    map2(x,y[names(x)], .f = function(a,b){
      a = a[b]
    })
  })

## Get networks from synapse
all.used.ids = c(all.used.ids, c('syn17015634', 'syn17015635'))
cor.mat = c(cohort1 = 'syn17015634', cohort2 = 'syn17015635') %>%
  map(.f =function(tblId){
    synapser::synTableQuery(paste0('select name,id from ', tblId))$asDataFrame() %>%
      dplyr::filter(name == 'Data_Correlation.txt') %>%
      dplyr::select(id) %>%
      dplyr::mutate(folderName = id) %>%
      mutate_at(.vars = 'folderName', .funs = function(id){
        map(id, .f = function(x){
          synGet(x, downloadFile = F)$properties$parentId
        }) %>% 
          map(.f = function(x){
            synGet(x, downloadFile = F)$properties$name
          }) %>% 
          unlist()
      }) %>%
      group_by(folderName) %>%
      nest() %>%
      deframe()
  })

## Get TF-DNA network from synapse
TF.DNA.net = fread(synGet('syn17060503')$path)
all.used.ids = c(all.used.ids, 'syn17060503')

bck.genes = map_dfr(dexp, function(x){
  map_dfr(x, function(y){y[,c('hgnc_symbol', 'ensembl_gene_id')]})
}) %>% unique()

TF.DNA.net = TF.DNA.net %>%
  tidyr::separate(from, c('from1'), sep = '\\_') %>%
  rename(from = from1) %>%
  dplyr::inner_join(bck.genes %>% rename(from = hgnc_symbol, from1 = ensembl_gene_id)) %>%
  dplyr::inner_join(bck.genes %>% rename(to = hgnc_symbol, to1 = ensembl_gene_id)) %>%
  dplyr::select(from1, to1) %>%
  distinct() %>%
  igraph::graph_from_data_frame(directed = T) %>%
  igraph::as_adjacency_matrix(type = 'both') %>%
  as.matrix()

## Run node ranking analysis
all.node.ranks = pmap(list(cor.mat, mod, dexp), .f = function(innerCmat, innerMod, innerDexp){
  tmp1 = map2(innerCmat, innerMod[names(innerCmat)], .f = function(stCmat, stMod){
    if(length(stMod) != 0){
      corMat = synGet(as.character(stCmat))$path %>%
        fread(data.table = FALSE, header = T) %>%
        dplyr::filter(row %in% unique(unlist(stMod)),
                      col %in% unique(unlist(stMod))) %>%
        igraph::graph_from_data_frame(directed = FALSE) %>%
        igraph::as_adjacency_matrix(attr = 'rho', type = 'both') %>%
        as.matrix()
      
      rnames = intersect(rownames(TF.DNA.net), rownames(corMat))
      cnames = intersect(colnames(TF.DNA.net), colnames(corMat))
      
      corMat[setdiff(rownames(corMat), rnames),
             setdiff(rownames(corMat), cnames)] = 0
      
      corMat[rnames, cnames] = corMat[rnames, cnames] * TF.DNA.net[rnames, cnames]
     
      nscores = future_map(stMod, .f = function(modGenes){
        corMat = corMat[modGenes, modGenes]
        g = igraph::graph_from_adjacency_matrix(corMat, mode = 'undirected',
                                                weighted = T, diag = F)
        scores = innerDexp$`Dx.Cranio+RNADaysToHarvest.Cranio-Control` %>%
          dplyr::full_join(variants$all %>%
                             dplyr::rename(hgnc_symbol = GeneName)) %>%
          dplyr::filter(ensembl_gene_id %in% modGenes)
        scores[is.na(scores)] = 0
        scores = scores %>%
          mutate(rk = rank(dscore) + rank(GERP) + rank(CADD) + rank(Polyphen2HDIV) + rank(Polyphen2HVAR),
                 rk = rank(rk)/max(rk)) %>%
          column_to_rownames('ensembl_gene_id') %>%
          as.data.frame()
        V(g)$Gx = scores[V(g)$name, 'rk']
        
        ns = scoreVertices(g, order = 3) 
        
        # Perform randomisation
        ns.perm = future_map(1:100, .f= function(i){
          V(g)$Gx = sample(scores$rk)
        }) %>%
          bind_cols()
        
        pval = map(1:length(ns), .f = function(i){
          sum(as.numeric(ns.perm[i,] >= ns[i]), na.rm = T)/100
        }) %>% 
          unlist() %>%
          p.adjust(method = 'BH')
        
        ns = data.frame(ensembl_gene_id = names(ns), nodeScores = ns, pval = pval) %>%
          left_join(scores %>% 
                      rownames_to_column(var = 'ensembl_gene_id')) %>%
          arrange(desc(nodeScores)) %>%
          dplyr::select(-rk)  
      })
    }
  })
}) %>%
  map(.f = function(x){
    x %>%
      map(.f = function(y){
        y %>% bind_rows(.id = 'moduleName')
      }) %>%
      bind_rows(.id = 'subType')
  }) %>%
  bind_rows(.id = 'cohort')

any.node.ranks = pmap(list(cor.mat, mod, dexp), .f = function(innerCmat, innerMod, innerDexp){
  tmp1 = map2(innerCmat, innerMod[names(innerCmat)], .f = function(stCmat, stMod){
    if(length(stMod) != 0){
      corMat = synGet(as.character(stCmat))$path %>%
        fread(data.table = FALSE, header = T) %>%
        dplyr::filter(row %in% unique(unlist(stMod)),
                      col %in% unique(unlist(stMod))) %>%
        igraph::graph_from_data_frame(directed = FALSE) %>%
        igraph::as_adjacency_matrix(attr = 'rho', type = 'both') %>%
        as.matrix()
      
      rnames = intersect(rownames(TF.DNA.net), rownames(corMat))
      cnames = intersect(colnames(TF.DNA.net), colnames(corMat))
      
      corMat[setdiff(rownames(corMat), rnames),
             setdiff(rownames(corMat), cnames)] = 0
      
      corMat[rnames, cnames] = corMat[rnames, cnames] * TF.DNA.net[rnames, cnames]
      
      nscores = future_map(stMod, .f = function(modGenes){
        corMat = corMat[modGenes, modGenes]
        g = igraph::graph_from_adjacency_matrix(corMat, mode = 'undirected',
                                                weighted = T, diag = F)
        scores = innerDexp$`Dx.Cranio+RNADaysToHarvest.Cranio-Control` %>%
          dplyr::full_join(variants$any %>%
                             dplyr::rename(hgnc_symbol = GeneName)) %>%
          dplyr::filter(ensembl_gene_id %in% modGenes)
        scores[is.na(scores)] = 0
        scores = scores %>%
          mutate(rk = rank(dscore) + rank(GERP) + rank(CADD) + rank(Polyphen2HDIV) + rank(Polyphen2HVAR),
                 rk = rank(rk)/max(rk)) %>%
          column_to_rownames('ensembl_gene_id') %>%
          as.data.frame()
        V(g)$Gx = scores[V(g)$name, 'rk']
        
        ns = scoreVertices(g, order = 3) 
        
        # Perform randomisation
        ns.perm = future_map(1:100, .f= function(i){
          V(g)$Gx = sample(scores$rk)
        }) %>%
          bind_cols()
        
        pval = map(1:length(ns), .f = function(i){
          sum(as.numeric(ns.perm[i,] >= ns[i]), na.rm = T)/100
        }) %>% 
          unlist() %>%
          p.adjust(method = 'BH')
        
        ns = data.frame(ensembl_gene_id = names(ns), nodeScores = ns, pval = pval) %>%
          left_join(scores %>% 
                      rownames_to_column(var = 'ensembl_gene_id')) %>%
          arrange(desc(nodeScores)) %>%
          dplyr::select(-rk)  
      })
    }
  })
}) %>%
  map(.f = function(x){
    x %>%
      map(.f = function(y){
        y %>% bind_rows(.id = 'moduleName')
      }) %>%
      bind_rows(.id = 'subType')
  }) %>%
  bind_rows(.id = 'cohort')
############################################################################################################

TF = fread(synGet('syn6040938')$path)
all.used.ids = c(all.used.ids,'syn6040938')

node.ranks = list(all = all.node.ranks, any = any.node.ranks) %>%
  bind_rows(.id = 'varianceComputation')
node.ranks$geneType = 'gene'
node.ranks$geneType[node.ranks$hgnc_symbol %in% TF$feature] = 'TF'

write.table(node.ranks, file = 'noderankingsTFDNA.tsv', sep = '\t')
obj = File('noderankingsTFDNA.tsv', name = 'Node Rankings (TF-DNA)', parentId = 'syn11635115')
obj = synStore(obj, executed = thisFile, use = all.used.ids, 
               activityName = activityName,
               activityDescription = activityDescription)