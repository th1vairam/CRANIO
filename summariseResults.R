# Get network
net.id = 'syn7342818'
all.used.ids = net.id
load(synGet(net.id)@filePath)

g = igraph::graph_from_adjacency_matrix(bicNetworks$network, mode = 'upper', diag = F)

# Query biomart to get gene names
mart <- useMart(host = "jul2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl", biomart = "ENSEMBL_MART_ENSEMBL")
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = V(g)$name, mart = mart)

# Get modules
mod = downloadFile('syn7347823') %>%
  dplyr::rename(ensembl_gene_id = Gene.ID) %>%
  left_join(Ensemble2HGNC)

mod.summary = mod %>%
  group_by(moduleLabel) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

tmp = filter(mod.summary, moduleLabel != 'turquoise')
p = ggplot(tmp, aes(count)) + geom_histogram()
p = p + theme_bw() + xlab('Module Size') + ylab('Number of Modules')
p

# Get enrichment results for modules
enrich.results = downloadFile('syn7460634') %>%
  filter(!is.na(Category))
ind1 = grep('Mus', enrich.results$GeneSetName)
ind2 = grep('mm9', enrich.results$GeneSetName)
ind3 = grep('mouse', enrich.results$GeneSetName)
ind4 = grep('MOUSE', enrich.results$GeneSetName)
ind = setdiff(1:dim(enrich.results)[1], c(ind1, ind2, ind3, ind4))
enrich.results = enrich.results[ind,]

tmp = enrich.results  %>%
  filter(Category == 'Cranio', fdr <= 0.05)
tmp = enrich.results  %>%
  filter(Category == 'CellTypeMarkers', fdr <= 0.05)
tmp = enrich.results  %>%
  filter(!(Category %in% c('CellTypeMarkers', 'Cranio')), fdr <= 0.05) %>% 
  ddply(.(ModuleLabel), .fun = function(x){
    x %>% 
      arrange(desc(Odds.Ratio)) %>% 
      top_n(5,Odds.Ratio)
    })

# Get cranio related gene list
cranio.gl = synGet('syn5752718')
load(cranio.gl@filePath)
cranio.gs = lapply(names(GeneSets), function(x, geneSets){
  as.data.frame(geneSets[[x]]) %>% 
    dplyr::mutate(y = TRUE) %>% 
    plyr::rename(c('geneSets[[x]]' = 'hgnc_symbol', 'y' = x))
}, GeneSets) %>%
  join_all(type = 'full')

# Get differential expression enriched gene scores
de.gene.scores = downloadFile('syn7477203') %>%
  dplyr::rename(ensembl_gene_id = Gene.ID) %>%
  left_join(mod) %>%
  left_join(cranio.gs) %>%
  dplyr::arrange(desc(global.regulator), desc(Scores))

# # Get variants enriched gene scores
# var.gene.scores = downloadFile('syn7477205') %>%
#   dplyr::rename(ensembl_gene_id = Gene.ID) %>%
#   left_join(mod) %>%
#   left_join(cranio.gs) %>%
#   dplyr::arrange(desc(global.regulator), desc(Scores))

# Get node degree
node.degree = igraph::degree(g) %>%
  rownameToFirstColumn('ensembl_gene_id') %>%
  dplyr::rename(node.degree=DF) %>%
  left_join(mod) %>%
  arrange(desc(node.degree))