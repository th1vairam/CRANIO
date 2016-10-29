#!/usr/bin/env Rscript

# Function to perform enrichment analysis of modules (from synapse as RData file)
# Get arguments from command line
args = c('syn7347823')

# Clear R console screen output
cat("\014")
############################################################################################################

############################################################################################################
#### Libraries ####
library(synapseClient)
library(CovariateAnalysis)
library(plyr)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(data.table)
library(biomaRt)

library(xlsx)

# Needs the dev branch
library(githubr)

# login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'enrichModules.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", 
                    ref="branch", 
                    refName='moduleEnrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0(thisFileName))

# Synapse specific parameters
activityName = 'Module enrichment'
activityDescription = 'Enrichment analysis of network modules using hypergeometric method'
############################################################################################################

############################################################################################################
# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                             genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genesInBackground # Background genes that are 
){
  genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
  genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
  genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
  genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
  
  pval = fisher.test(
    matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
             length(intersect(genesInGeneSet, genesInNonSignificantSet)),
             length(intersect(genesOutGeneSet, genesInSignificantSet)),
             length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
           nrow=2, ncol=2),
    alternative="greater")
  OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genesInGeneSet),
                    noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                    Odds.Ratio = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')
  )
  )
}
############################################################################################################

############################################################################################################
#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath)

gsets = c("GO_Biological_Process", "KEGG_2015", "Reactome_2015", "WikiPathways_2015")
GeneSets.Enrichr = GeneSets[gsets]

# Download AD related gene sets from synapse
GL_OBJ = synGet('syn4893059');
ALL_USED_IDs = c(ALL_USED_IDs, GL_OBJ$properties$id)
load(GL_OBJ@filePath)

GeneSets.CM = list(CellTypeMarkers = GeneSets[grep('Zhang', names(GeneSets))])

# Download cranio related gene sets from synapse
GL_OBJ = synGet('syn5752718');
ALL_USED_IDs = c(ALL_USED_IDs, GL_OBJ$properties$id)
load(GL_OBJ@filePath)

GeneSets.Cranio = list(Cranio = GeneSets)
############################################################################################################

############################################################################################################
#### Get modules ####
# Download modules from synapse
MOD_OBJ = synapseClient::synGet(args[1])
ALL_USED_IDs = c(ALL_USED_IDs, MOD_OBJ$properties$id)
FNAME = str_replace_all(tools::file_path_sans_ext(MOD_OBJ$properties$name), 'Modules', 'Enrichment')
parentId = MOD_OBJ$properties$parentId

# Load modules 
MOD = data.table::fread(MOD_OBJ@filePath, data.table=F, header=T)
############################################################################################################

############################################################################################################
#### Background gene list ####
# Convert ensemble gene id's to hgnc symbols using biomart
ensembl = useMart(host = "jul2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl", biomart = "ENSEMBL_MART_ENSEMBL")
ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$Gene.ID, mart = ensembl)
backGroundGenes = unique(ensg2hgnc$hgnc_symbol)

MOD = merge(MOD, ensg2hgnc, by.x = 'Gene.ID', by.y = 'ensembl_gene_id', all.x=T)
############################################################################################################

############################################################################################################
#### Filter gene list ####
GeneSets = c(GeneSets.Enrichr, GeneSets.CM, GeneSets.Cranio)
GeneSets = filterGeneSets(GeneSets, backGroundGenes, minSize = 10, maxSize = 1000)
############################################################################################################

############################################################################################################
#### Perform enrichment analysis ####
# Perform enrichment analysis (for modules greater than 20 genes only)
enrichResults = list()
for (name in unique(MOD$moduleLabel)){
  genesInModule = unique(MOD$hgnc_symbol[MOD$moduleLabel == name])  
  if (length(genesInModule) > 20){
    enrichResults[[name]] = ldply(GeneSets,
                                   .fun = function(x, genesInModule, genesInBackground, fisherEnrichment){
                                     ldply(x, .fun = fisherEnrichment, genesInModule, genesInBackground, .parallel = T, .id = 'GeneSetName')
                                     }, 
                                  unique(genesInModule), unique(backGroundGenes), fisherEnrichment,
                                  .parallel = T, .id = 'Category') %>%
      dplyr::mutate(fdr = p.adjust(pval, 'fdr'))
    } else {
      enrichResults[[name]] = data.frame(Category = NA, GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, Odds.Ratio = NA, Genes = NA, fdr = NA)
  }
  writeLines(paste0('Completed ',name))  
}

# Write results to file
tmp = rbindlist(enrichResults, use.names = TRUE, idcol = 'ModuleLabel', fill = TRUE) %>%
  filter(!is.na(Category))

write.table(tmp, file = paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), sep='\t', row.names=F)
gc()
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
ENR_OBJ = File(paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), name = paste(FNAME,'FishersExact'), parentId = parentId)
annotations(ENR_OBJ) = list(fileType = 'tsv',
                            algorithm = 'FishersExact',
                            resultsType = 'moduleEnrichment')
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Filter modules with more than 20 genes
filteredModules = MOD %>%
  group_by(moduleLabel) %>%
  dplyr::summarise(count = length(unique(hgnc_symbol))) %>%
  filter(count >= 20)

# Re-assign smaller modules to 0/NoModule
ind = which(!(MOD$moduleLabel %in% filteredModules$moduleLabel))
MOD[ind, 'moduleNumber'] = 0
MOD[ind, 'moduleLabel'] = 'NoModule'
############################################################################################################

############################################################################################################
# Filter enrichment results per module and present the pathways
# Remove mouse related gene sets
enrichResults = rbindlist(enrichResults, use.names = TRUE, idcol = 'ModuleLabel', fill = TRUE) %>%
  filter(!is.na(Category))
ind1 = grep('Mus', enrichResults$GeneSetName)
ind2 = grep('mm9', enrichResults$GeneSetName)
ind3 = grep('mouse', enrichResults$GeneSetName)
ind4 = grep('MOUSE', enrichResults$GeneSetName)
ind = setdiff(1:dim(enrichResults)[1], c(ind1, ind2, ind3, ind4))
enrichResults = enrichResults[ind,]

# Filter enrichment results (pathways)
tmp = dlply(enrichResults, .(ModuleLabel),.fun = function(x){
  x1 = x %>%
    filter(Category %in% c("CellTypeMarkers", "Cranio")) %>%
    mutate(fdr = p.adjust(pval, 'fdr')) %>%
    filter(fdr <= 0.05) %>% arrange(desc(Odds.Ratio))
  return(x1)
})
tmp = tmp[sapply(tmp, dim)[1,] != 0]

for(mod in names(tmp))
  write.xlsx(tmp[[mod]], 
             file = 'cranioPathwaysCaseControl.xlsx', 
             sheetName = paste(mod,'Pathways'), 
             row.names=F, col.names=T, append = T)

# Write results to synapse
OBJ = File('cranioPathwaysCaseControl.xlsx', name = 'Cranio Enrichment Pathways Filtered (case+control)', parentId = 'syn7320952')
OBJ = synStore(OBJ, executed = thisFile, used = ENR_OBJ, activityName = activityName,
               activityDescription = activityDescription)
############################################################################################################