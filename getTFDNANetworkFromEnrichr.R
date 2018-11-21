############################################################################################################
#### Curate TF-DNA network from Enrichr ####
net.urls = c('http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
             'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_TF_ChIP-seq_2015',
             'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=TF_Perturbations_Followed_by_Expression')

encode.gs = map_dfr(net.urls, .f = function(urlId){
  encode.data = url(urlId) %>%
    read.csv()
  
  gs = map_dfr(1:dim(encode.data)[1], .f = function(x, tmp){
    genes = str_split(tmp[x,1],'\t')[[1]]
    data.frame(to = genes[3:length(genes)]) %>%
      mutate(from = genes[1]) 
  }, encode.data) %>%
    dplyr::select(from,to) %>%
    distinct()
}) %>%
  distinct()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'getTFDNANetworkFromEnrichr.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CRANIO", ref="branch", refName='rankNodes')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0(thisFileName))

# Synapse specific parameters
activityName = 'Enrichr gene set curation'
############################################################################################################

############################################################################################################
#### Store in synapse ####
write.table(encode.gs, file = 'encodeNet.tsv', sep = '\t', row.names = F)
obj = File('encodeNet.tsv', name = 'Enrichr TF-DNA Network', parentId = 'syn11635115')
obj = synStore(obj, executed = thisFile, used = net.urls, activityName = activityName)

############################################################################################################