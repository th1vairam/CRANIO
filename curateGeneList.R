############################################################################################################
# Function to curate important genelists for enrichment

# Clear R console screen output
cat("\014")

# Libraries 
library(synapseClient)
library(githubr)
library(CovariateAnalysis)
library(dplyr)
library(data.table)
library(biomaRt)
library(tidyr)

# Login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'curateGeneList.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", 
                    ref="branch", 
                    refName='moduleEnrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0(thisFileName))

# Synapse specific parameters
activityName = 'Geneset curation'
activityDescription = 'Curating CRANIO related important genesets'
############################################################################################################

############################################################################################################
#### Differentially expressed genes ####
diffexp = downloadFile('syn5745241') %>%
  dplyr::rename(ensembl_gene_id = V1)

# Convert ensemble gene id's to hgnc symbols using biomart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
                  values = diffexp$ensembl_gene_id, mart = ensembl)

# Get gene hgnc_symbols
diffexp = left_join(diffexp, ensg2hgnc)

GeneSets = list()
GeneSets[['Case_vs_Control_0.05FDR_2FC']] = diffexp %>%
  filter(adj.P.Val <= 0.05, abs(logFC) >= 1) %>%
  dplyr::select(hgnc_symbol) %>% unique %>% unlist %>%
  as.character

#### Differentially expressed genes in clusters ####
load(synGet('syn4939237')@filePath)
GeneSets = c(GeneSets, 
             lapply(DEgeneLists, function(x, diffExp){
               diffExp %>%
                 filter(ensembl_gene_id %in% x) %>%
                 dplyr::select(hgnc_symbol) %>% unique %>% unlist %>%
                 as.character
              }, diffexp))

#### Differentially expressed genes in clusters ####
impgl = downloadFile('syn3105988') %>%
  gather(Pathway, Presence, -V1) %>%
  group_by(Pathway, Presence) %>%
  summarise(V1 = list(unique(sort(V1)))) %>%
  filter(Presence == 'YES')
tmp = impgl$V1
names(tmp) = impgl$Pathway

GeneSets = c(GeneSets, tmp)
############################################################################################################

############################################################################################################
#### Write gene list to synapse ####
# Write results to synapse
save(list = 'GeneSets', file = 'CranioGeneSets.RData')

OBJ = File('CranioGeneSets.RData', name = 'Cranio Gene Sets', parentId = 'syn5650456')
OBJ = synStore(OBJ, executed = thisFile, used = c('syn3105988', 'syn4939237', 'syn5745241'),
               activityName = activityName,
               activityDescription = activityDescription)
############################################################################################################