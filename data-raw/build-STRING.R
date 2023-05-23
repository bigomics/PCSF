
getSTRING <- function(min.score=800) {
  require(STRINGdb)
  require(igraph)
  require(biomaRt)

  # 1. getSTRINGdb for human
  string_db <- STRINGdb$new(species=9606, score_threshold=min.score)
  graph <- string_db$get_graph()
  graph
    
  # 4. map gene ids to protein ids
  ### get gene/protein ids via Biomart
  mart=useMart(##host = 'https://grch37.ensembl.org',
    biomart='ENSEMBL_MART_ENSEMBL',
    dataset='hsapiens_gene_ensembl')
  
  ### extract protein ids from the human network
  protein_ids <- sapply(strsplit(V(graph)$name, '\\.'),
    function(x) x[2])
  
  ### get protein to gene id mappings
  mart_results <- getBM(attributes = c("hgnc_symbol",
    "ensembl_peptide_id"),
    filters = "ensembl_peptide_id", values = protein_ids,
    mart = mart)
  head(mart_results)
  
  ### replace protein ids with gene ids
  ix <- match(protein_ids, mart_results$ensembl_peptide_id)
  gene <- mart_results$hgnc_symbol[ix]
  
  V(graph)$name <- gene
  has.name <- !V(graph)$name %in% c(NA,"")
  graph <- subgraph(graph, which(has.name))
  graph
  
  ee <- get.edgelist(graph)
  STRING <- data.frame( ee, 1000-E(graph)$combined_score)
  colnames(STRING) <- c("from","to","cost")
  STRING <- STRING[which(!is.na(STRING$from) & !is.na(STRING$to)),]
  STRING$cost <- STRING$cost / max(STRING$cost,na.rm=TRUE)
  head(STRING)
  STRING
}

STRING <- getSTRING(750)
dim(STRING)
usethis::use_data(STRING, overwrite = TRUE)
