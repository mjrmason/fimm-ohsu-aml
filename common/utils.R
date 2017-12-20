invert.list <- function(lst) {
  tbl <- ldply(1:length(lst), .parallel = FALSE,
               .fun = function(i) {
                  if(is.null(lst[[i]])) { return(NULL) }
                  df <- data.frame(gene = lst[[i]])
                  df$gene.set <- names(lst)[i]
                  df
               })
  dlply(tbl, .variables = "gene", .fun = function(df) df$gene.set)
}

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

ensg.to.sym.mapping <- function(gene.ids) {
  # Get a mapping from ensembl id to hugo symbol
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'ensembl_gene_id', 
              values = gene.ids,
              mart = ensembl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$symbol %in% c("")),]
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

symbols.to.ensg.mapping <- function(symbols, ensembl = NULL, curl = NULL) {
  # Get a mapping from gene symbol to ensembl id
  if(is.null(ensembl)) {
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  }
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl,
              curl = curl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}
