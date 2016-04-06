#!/usr/bin/env Rscript

# Function to perform enrichment analysis of modules (from synapse as RData file)
# Get arguments from command line
args = c('syn5700963')

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

gsets = c("Achilles_fitness_decrease", "Achilles_fitness_increase", "Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up",
          "BioCarta_2015", "CMAP_down", "CMAP_up", "Chromosome_Location", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Cross_Species_Phenotype",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2015", "ESCAPE", "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function",
          "GeneSigDB", "HMDB_Metabolites", "HomoloGene", "Human_Gene_Atlas", "KEGG_2015",
          "MGI_Mammalian_Phenotype", "MGI_Mammalian_Phenotype_Level_3", "MGI_Mammalian_Phenotype_Level_4",
          "MSigDB_Oncogenic_Signatures", "Mouse_Gene_Atlas", "NCI-60_Cancer_Cell_Lines", "NCI-Nature", 
          "NURSA_Human_Endogenous_Complexome", "OMIM_Disease", "OMIM_Expanded", "PPI_Hub_Proteins", "Pfam_InterPro_Domains",
          "Phosphatase_Substrates_from_DEPOD", "Reactome_2015", "SILAC_Phosphoproteomics", "TF-LOF_Expression_from_GEO", 
          "TargetScan_microRNA", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
          "Transcription_Factor_PPIs", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "WikiPathways_2015",
          "GeneFamily","CellMarkers")
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
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$GeneIDs, mart = ensembl)
backGroundGenes = unique(ensg2hgnc$hgnc_symbol)

MOD = merge(MOD, ensg2hgnc, by.x = 'GeneIDs', by.y = 'ensembl_gene_id', all.x=T)
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
for (name in unique(MOD$modulelabels)){
  genesInModule = unique(MOD$hgnc_symbol[MOD$modulelabels == name])  
  if (length(genesInModule) > 20){
    enrichResults[[name]] = lapply(GeneSets,
                                   function(x, genesInModule, genesInBackground){
                                     tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                                     tmp = rownameToFirstColumn(tmp,'GeneSetName')
                                     return(tmp)
                                     },
                                   unique(genesInModule), unique(backGroundGenes)) %>%
      rbindlist(use.names=TRUE, idcol = 'Category') %>%
      dplyr::mutate(fdr = p.adjust(pval, 'fdr'))
  } else {
    enrichResults[[name]] = data.frame(GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, OR = NA, Category = NA, fdr = NA)
  }
  writeLines(paste0('Completed ',name))  
}

# Write results to file
tmp = rbindlist(enrichResults, use.names = TRUE, idcol = 'ComparisonName', fill = TRUE) %>%
  filter(!is.na(Category))

tmp$pval = unlist(tmp$pval)
tmp$ngenes = unlist(tmp$ngenes)
tmp$noverlap = unlist(tmp$noverlap)
tmp$Odds.Ratio = unlist(tmp$Odds.Ratio)
tmp$Genes = unlist(tmp$Genes)

write.table(tmp, file = paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), sep='\t', row.names=F)
gc()
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
algo = 'Fisher'
ENR_OBJ = File(paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), name = paste(FNAME,algo), parentId = parentId)
annotations(ENR_OBJ) = annotations(MOD_OBJ)
ENR_OBJ@annotations$fileType = 'tsv'
ENR_OBJ@annotations$enrichmentMethod = 'Fisher'
ENR_OBJ@annotations$enrichmentGeneSet = 'Enrichr and AD'
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Filter modules with more than 20 genes
filteredModules = MOD %>%
  group_by(modulelabels) %>%
  dplyr::summarise(count = length(unique(hgnc_symbol))) %>%
  filter(count >= 20)

# Re-assign smaller modules to 0/NoModule
ind = which(!(MOD$modulelabels %in% filteredModules$modulelabels))
MOD[ind, 'moduleNumber'] = 0
MOD[ind, 'modulelabels'] = 'NoModule'
############################################################################################################

############################################################################################################
# Calculate node degree within module as well as in the whole network
# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet('syn5700531')
ALL_USED_IDs = c(ALL_USED_IDs, 'syn5700531')

# Load sparse network
load(NET_OBJ@filePath)

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(sparseNetwork, mode = 'undirected', weighted = T, diag = F)

# Calculate node degree and merge with modules
MOD = MOD %>%
  left_join( degree(g, v = V(g), mode = "all") %>%
               rownameToFirstColumn('GeneIDs')) %>%
  dplyr::rename(node.degree = DF)

# Calculate module specific node degree
tmp = MOD %>% 
  dplyr::select(GeneIDs, modulelabels, hgnc_symbol, node.degree) %>%
  plyr::dlply(.(modulelabels), .fun = function(x){
    sg = subgraph(g, x$GeneIDs)
    x = left_join(x, 
                  degree(sg, v = V(sg), mode = "all") %>%
                    rownameToFirstColumn('GeneIDs')) %>%
      dplyr::rename(node.degree.in.module = DF) %>%
      arrange(desc(node.degree.in.module))
})

# Store node degree in file
for(mod in names(tmp))
  write.xlsx(tmp[[mod]], 
             file = 'cranioNodeDegree.xlsx', 
             sheetName = paste(mod,'Degree'), row.names=F, col.names=F, append = T)

# Write results to synapse
OBJ = File('cranioNodeDegree.xlsx', name = 'Cranio Node Degree', parentId = 'syn5879683')
OBJ = synStore(OBJ, executed = thisFile, used = c('syn5700531', args[1]), activityName = activityName,
               activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Filter enrichment results per module and present the pathways
enrichResults = downloadFile('syn5701314')

# Remove mouse related gene sets
ind1 = grep('Mus', enrichResults$GeneSetName)
ind2 = grep('mm9', enrichResults$GeneSetName)
ind3 = grep('mouse', enrichResults$GeneSetName)
ind4 = grep('MOUSE', enrichResults$GeneSetName)
ind = setdiff(1:dim(enrichResults)[1], c(ind1, ind2, ind3, ind4))
enrichResults = enrichResults[ind,]

# Filter enrichment results (pathways)
tmp = dlply(enrichResults, .(ComparisonName),.fun = function(x){
  x1 = x %>%
    filter(Category %in% c("WikiPathways_2015", "CellTypeMarkers", "Cranio", "Reactome_2015", "GO_Biological_Process")) %>%
    mutate(fdr = p.adjust(pval, 'fdr')) %>%
    filter(fdr <= 0.05) %>% arrange(desc(Odds.Ratio))
  x1$OR = NULL
  return(x1)
})
tmp = tmp[sapply(tmp, dim)[1,] != 0]

for(mod in names(tmp))
  write.xlsx(tmp[[mod]], 
             file = 'cranioPathways.xlsx', 
             sheetName = paste(mod,'Pathways'), 
             row.names=F, col.names=T, append = T)

# Write results to synapse
OBJ = File('cranioPathways.xlsx', name = 'Cranio Pathways', parentId = 'syn5879683')
OBJ = synStore(OBJ, executed = thisFile, used = c('syn5701314'), activityName = activityName,
               activityDescription = activityDescription)

# Filter enrichment results (regulations)
tmp = dlply(enrichResults, .(ComparisonName),.fun = function(x){
  x1 = x %>%
    filter(Category %in% c("TargetScan_microRNA", "TF-LOF_Expression_from_GEO", "ChEA", "Chromosome_Location")) %>%
    mutate(fdr = p.adjust(pval, 'fdr')) %>%
    filter(fdr <= 0.05) %>% arrange(desc(Odds.Ratio))
  x1$OR = NULL
  return(x1)
})
tmp = tmp[sapply(tmp, dim)[1,] != 0]

for(mod in names(tmp))
  write.xlsx(tmp[[mod]], 
             file = 'cranioRegulations.xlsx', 
             sheetName = paste(mod,'TFs miRs Chromosome'), 
             row.names=F, col.names=T, append = T)

# Write results to synapse
OBJ = File('cranioRegulations.xlsx', name = 'Cranio TFs miRNAs Chromosomes', parentId = 'syn5879683')
OBJ = synStore(OBJ, executed = thisFile, used = c('syn5701314'), activityName = activityName,
               activityDescription = activityDescription)
############################################################################################################