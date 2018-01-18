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
##  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "may2017.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'ensembl_gene_id', 
              values = gene.ids,
              mart = ensembl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$symbol %in% c("")),]
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

symbols.to.ensg.mapping <- function(symbols) {
  dummy <- data.frame(query = symbols, ensembl = NA)
  bm <- tryCatch({queryMany(symbols, scopes="symbol", fields=c("ensembl.gene"), species="human")}, error = function(e) { return(dummy) })
  flag <- grepl(pattern="ensembl", colnames(bm))
  if(length(which(flag)) != 1) {
    stop(paste0("Could not find ensembl col in: ", paste(colnames(bm), collapse=" "), "\n"))
  }
  ensg.col <- colnames(bm)[flag]
  bm <- bm[, c("query", ensg.col)]
  lst <- bm[,ensg.col]
  names(lst) <- bm$query
  bm <- ldply(lst, .fun = function(comp) data.frame(unlist(comp)))
  colnames(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),,drop=F]
  bm <- bm[!is.na(bm$ensg) & !is.na(bm$symbol),,drop=F]
  bm
}

symbols.to.ensg.mapping.biomart <- function(symbols, ensembl = NULL, curl = NULL) {
  # Get a mapping from gene symbol to ensembl id
  if(is.null(ensembl)) {
##    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "may2017.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
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
